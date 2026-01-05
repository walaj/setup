# app.R
# Minimal Shiny console to step through loci, drive desktop IGV, and log decisions.

library(shiny)
library(data.table)
library(httr)

# --- CONFIG -------------------------------------------------------------
IGV_HOST <- "http://127.0.0.1"
IGV_PORT <- 60151 # Desktop IGV: View -> Preferences -> Advanced -> Enable port (60151)
BASE <- sprintf("%s:%d", IGV_HOST, IGV_PORT)

library(httr)

igv_cmd <- function(path, base = sprintf("http://127.0.0.1:%d", 60151), timeout_sec = 10) {
  url <- paste0(base, "/", path)
  r <- try(GET(url, timeout(timeout_sec)), silent = TRUE)
  if (inherits(r, "try-error")) return(list(ok = FALSE, text = "ERROR timeout/connection"))
  txt <- try(content(r, "text", encoding = "UTF-8"), silent = TRUE)
  if (inherits(txt, "try-error")) txt <- ""
  list(ok = !grepl("^ERROR", txt, ignore.case = TRUE), text = txt)
}

igv_goto <- function(locus_string, base = sprintf("http://127.0.0.1:%d", 60151)) {
  res <- igv_cmd(paste0("goto?locus=", URLencode(locus_string)), base = base)
  res$ok
}

# Accepts several param names (IGV versions differ)
igv_set_snapshot_dir <- function(dir, base = sprintf("http://127.0.0.1:%d", 60151)) {
  dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
  for (parm in c("directory","dir","path")) {
    if (igv_cmd(paste0("snapshotDirectory?", parm, "=", URLencode(dir)), base = base)$ok) return(TRUE)
  }
  FALSE
}

# Try both 'filename' and 'file'
igv_snapshot <- function(filename, base = sprintf("http://127.0.0.1:%d", 60151)) {
  for (parm in c("filename","file")) {
    if (igv_cmd(paste0("snapshot?", parm, "=", URLencode(filename)), base = base)$ok) return(TRUE)
  }
  FALSE
}

igv_get_locus <- function(...) NA_character_

# Try to query current locus from IGV (if supported by your IGV build).
# If your IGV doesn’t expose this, we’ll fall back to the expected locus.
# igv_get_locus <- function() {
#   url <- paste0(BASE, "/locus")
#   out <- try(httr::content(httr::GET(url, timeout(2)), as = "text"), silent = TRUE)
#   if (inherits(out, "try-error")) return(NA_character_)
#   # Some IGV builds return plain text; others return JSON. Keep it simple:
#   gsub("\\s+", "", out)
# }

# --- APP ---------------------------------------------------------------
ui <- fluidPage(
  tags$head(
    # Minimal hotkeys (j/k or left/right arrows)
    tags$script(HTML("
document.addEventListener('keydown', function(e) {
if (e.key === 'j' || e.key === 'ArrowRight') Shiny.setInputValue('hotkey_next', Date.now());
if (e.key === 'k' || e.key === 'ArrowLeft') Shiny.setInputValue('hotkey_prev', Date.now());
if (e.key === '1') Shiny.setInputValue('hotkey_good', Date.now());
if (e.key === '2') Shiny.setInputValue('hotkey_borderline', Date.now());
if (e.key === '3') Shiny.setInputValue('hotkey_bad', Date.now());
if (e.key === 's') Shiny.setInputValue('hotkey_snapshot', Date.now());
});
"))
  ),
  titlePanel("IGV Review Console"),
  sidebarLayout(
    sidebarPanel(
      fileInput("loci_file", "Load loci (TSV/CSV with chr,start,end[,name])", accept = c(".tsv",".csv",".txt")),
      textInput("genome", "Genome label (for your log only)", "hg38"),
      numericInput("port", "IGV port", IGV_PORT, min = 1, max = 65535, step = 1),
      textInput("snap_dir", "Snapshot directory (optional)", ""),
      textInput("session_id", "Session label (used in filenames/log)", format(Sys.time(), "%Y%m%d_%H%M%S")),
      actionButton("connect", "Connect / Initialize", class = "btn-primary"),
      hr(),
      h4("Controls"),
      actionButton("prev", "◀︎ Prev (k / ←)"),
      actionButton("next_btn", "Next ▶︎ (j / →)", class = "btn-success"),
      br(), br(),
      radioButtons("call", "Your call",
                   choices = c("Good [1]" = "good", "Borderline [2]" = "borderline", "Bad [3]" = "bad"),
                   selected = character(0), inline = TRUE),
      textInput("note", "Optional note", ""),
      actionButton("save", "Save call + Next", class = "btn-primary"),
      actionButton("snapshot", "Snapshot (s)"),
      br(), br(),
      actionButton("backfix", "Save call only (don’t advance)"),
      hr(),
      checkboxInput("auto_snapshot_on_save", "Auto snapshot on Save", TRUE),
      helpText("Hotkeys: j/→=Next, k/←=Prev, 1/2/3=set call, s=Snapshot")
    ),
    mainPanel(
      h4(textOutput("status")),
      verbatimTextOutput("current_locus_txt"),
      tableOutput("progress_table")
    )
  )
)
server <- function(input, output, session) {
  library(httr)
  
  # ---------------- IGV helpers (quiet + tolerant) ----------------
  base_url <- reactiveVal(sprintf("http://127.0.0.1:%d", IGV_PORT))
  observeEvent(input$port, ignoreInit = TRUE, {
    base_url(sprintf("http://127.0.0.1:%d", as.integer(input$port)))
  })
  
  igv_cmd <- function(path, timeout_sec = 10) {
    url <- paste0(base_url(), "/", path)
    r <- try(GET(url, timeout(timeout_sec)), silent = TRUE)
    if (inherits(r, "try-error")) return(list(ok = FALSE, text = "ERROR timeout/connection"))
    txt <- try(content(r, "text", encoding = "UTF-8"), silent = TRUE)
    if (inherits(txt, "try-error")) txt <- ""
    list(ok = (httr::status_code(r) == 200L && !grepl("^ERROR", txt, ignore.case = TRUE)),
         text = txt)
  }
  
  igv_goto <- function(locus_string) {
    res <- igv_cmd(paste0("goto?locus=", utils::URLencode(locus_string, reserved = TRUE)))
    res$ok
  }
  
  # Accept multiple param names across IGV versions
  igv_set_snapshot_dir <- function(dir) {
    dir <- normalizePath(dir, winslash = "/", mustWork = FALSE)
    for (parm in c("directory","dir","path")) {
      if (igv_cmd(paste0("snapshotDirectory?", parm, "=", utils::URLencode(dir, reserved = TRUE)))$ok)
        return(TRUE)
    }
    FALSE
  }
  
  # Try both 'filename' and 'file'
  igv_snapshot <- function(filename) {
    for (parm in c("filename","file")) {
      if (igv_cmd(paste0("snapshot?", parm, "=", utils::URLencode(filename, reserved = TRUE)))$ok)
        return(TRUE)
    }
    FALSE
  }
  
  # Many IGV builds don't expose /locus; return NA to keep logs clean
  igv_get_locus <- function() NA_character_
  
  # ---------------- Reactive state ----------------
  loci_dt   <- reactiveVal(NULL)
  idx       <- reactiveVal(1L)
  log_dt    <- reactiveVal(NULL)
  snap_ready <- reactiveVal(FALSE)
  
  # ---------------- Load loci file ----------------
  observeEvent(input$loci_file, {
    ext <- tools::file_ext(input$loci_file$name)
    dt <- tryCatch({
      if (tolower(ext) %in% "csv") data.table::fread(input$loci_file$datapath)
      else data.table::fread(input$loci_file$datapath, sep = "\t")
    }, error = function(e) NULL)
    validate(need(!is.null(dt), "Could not read loci file."))
    
    # Normalize columns
    cn <- tolower(names(dt))
    data.table::setnames(dt, cn)
    need_cols <- c("chr","start","end")
    validate(need(all(need_cols %in% names(dt)),
                  "File must have columns: chr, start, end (and optional name)"))
    
    if (!"name" %in% names(dt)) dt[, name := sprintf("%s:%s-%s", chr, start, end)]
    dt[, locus := sprintf("%s:%s-%s", chr, start, end)]
    dt[, id := .I]
    
    loci_dt(dt)
    idx(1L)
    log_dt(NULL)
    
    # Jump to first locus on load
    row <- dt[1L]
    if (nrow(dt) > 0L) {
      ok <- igv_goto(row$locus)
      if (!isTRUE(ok)) showNotification("IGV navigation failed (port enabled? correct port?).", type = "error")
    }
  })
  
  # ---------------- Connect / init ----------------
  observeEvent(input$connect, {
    if (nzchar(input$snap_dir)) {
      ok <- igv_set_snapshot_dir(input$snap_dir)
      snap_ready(ok)
      showNotification(if (ok) "Snapshot directory set in IGV."
                       else "Failed to set snapshot directory. Check path & IGV port.",
                       type = if (ok) "message" else "error")
    } else {
      snap_ready(FALSE)
      showNotification("No snapshot directory set. You can still review; snapshots disabled.", type = "warning")
    }
  })
  
  # ---------------- Helpers ----------------
  current_row <- reactive({
    dt <- loci_dt()
    if (is.null(dt) || nrow(dt) < 1) return(NULL)
    i <- max(1L, min(idx(), nrow(dt)))
    dt[i]
  })
  
  goto_current_in_igv <- function() {
    row <- current_row(); req(row)
    ok <- igv_goto(row$locus)
    if (!isTRUE(ok)) showNotification("IGV navigation failed (port enabled? correct port?).", type = "error")
  }
  
  move_idx <- function(delta) {
    dt <- loci_dt(); req(dt)
    new_i <- max(1L, min(idx() + delta, nrow(dt)))
    idx(new_i)
    goto_current_in_igv()
  }
  
  # ---------------- Navigation ----------------
  observeEvent(input$prev,          { move_idx(-1) })
  observeEvent(input[["next"]],     { move_idx(+1) })  # 'next' is a keyword; use [["next"]]
  observeEvent(input$hotkey_prev,   { move_idx(-1) })
  observeEvent(input$hotkey_next,   { move_idx(+1) })
  
  # Hotkeys to set calls
  observeEvent(input$hotkey_good,       { updateRadioButtons(session, "call", selected = "good") })
  observeEvent(input$hotkey_borderline, { updateRadioButtons(session, "call", selected = "borderline") })
  observeEvent(input$hotkey_bad,        { updateRadioButtons(session, "call", selected = "bad") })
  observeEvent(input$hotkey_snapshot,   { isolate({ if (snap_ready()) take_snapshot() }) })
  
  # ---------------- Snapshot wrapper ----------------
  take_snapshot <- function() {
    row <- current_row(); req(row)
    if (!snap_ready()) {
      showNotification("Snapshot dir not set. Click Connect and set a snapshot directory.", type = "warning")
      return(invisible(FALSE))
    }
    fname <- sprintf("%s_%s_%s_%06d.png",
                     input$session_id, input$genome, row$name, row$id)
    ok <- igv_snapshot(fname)
    if (ok) showNotification(paste("Snapshot saved:", fname)) else
      showNotification("Snapshot failed.", type = "error")
    invisible(ok)
  }
  observeEvent(input$snapshot, { take_snapshot() })
  
  # ---------------- Logging ----------------
  append_log <- function(call_val, note_val, locus_seen = NA_character_) {
    row <- current_row(); req(row)
    seen <- locus_seen
    if (is.na(seen) || !nzchar(seen)) seen <- igv_get_locus()
    if (is.na(seen) || !nzchar(seen)) seen <- row$locus
    
    dt <- data.table::data.table(
      session   = input$session_id,
      genome    = input$genome,
      id        = row$id,
      name      = row$name,
      chr       = row$chr,
      start     = row$start,
      end       = row$end,
      locus     = row$locus,   # planned locus
      igv_locus = seen,        # IGV-reported if available
      call      = call_val,
      note      = note_val,
      ts        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )
    cur <- log_dt()
    log_dt(if (is.null(cur)) dt else data.table::rbindlist(list(cur, dt), use.names = TRUE))
  }
  
  # Save + (optional) snapshot + advance
  observeEvent(input$save, {
    req(input$call)
    if (isTRUE(input$auto_snapshot_on_save) && snap_ready()) take_snapshot()
    append_log(input$call, input$note)
    
    # Clear selection for next row
    updateRadioButtons(session, "call", selected = character(0))
    updateTextInput(session, "note", value = "")
    
    # Persist CSV (atomic-ish)
    out <- file.path(getwd(), sprintf("igv_review_%s.csv", input$session_id))
    data.table::fwrite(log_dt(), out)
    
    move_idx(+1)
  })
  
  # Save without advancing
  observeEvent(input$backfix, {
    req(input$call)
    append_log(input$call, input$note)
    out <- file.path(getwd(), sprintf("igv_review_%s.csv", input$session_id))
    data.table::fwrite(log_dt(), out)
    showNotification("Saved without advancing.")
  })
  
  # ---------------- UI status ----------------
  output$status <- renderText({
    dt <- loci_dt()
    if (is.null(dt)) return("Load a loci file to begin.")
    paste0("Reviewing ", idx(), " / ", nrow(dt))
  })
  
  output$current_locus_txt <- renderText({
    row <- current_row()
    if (is.null(row)) return("")
    paste("Current locus:", row$locus, "| name:", row$name)
  })
  
  output$progress_table <- renderTable({
    ld <- log_dt()
    if (is.null(ld)) return(NULL)
    tail(ld[, .(id, name, call, igv_locus, ts, note)], 10)
  })
}

shinyApp(ui, server)
