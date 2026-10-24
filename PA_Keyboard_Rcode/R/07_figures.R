# =============================================================================
# Publication figures: scenario truths and dose-finding performance
# =============================================================================
# This module creates the two manuscript figures as TIFF files. Figure 1
# displays truth-level adult/pediatric exposure correspondence and toxicity
# compatibility; Figure 2 displays simulation-based PCS and overtoxic allocation.
#
# Function guide:
#   create_figure1()                     Plot scenario AUC/DLT relationships.
#   create_figure2()                     Plot dose-finding operating characteristics.
#   Plot-specific helper functions       Defined locally within each figure routine.
# =============================================================================

# =============================================================================
# Figure 1: scenario-level exposure matching and toxicity compatibility
# =============================================================================
create_figure1 <- function(output_tiff) {
  if (!requireNamespace("ggplot2", quietly  = TRUE)) {
    stop("Package 'ggplot2' is required to create Figure 1.", call.  = FALSE)
  }
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(dirname(output_tiff), recursive  = TRUE, showWarnings  = FALSE)
  # -----------------------------------------------------------------------------

  AUC_CALIPER <- 0.22
  AUC_TIE_TOL <- 1e-12
  DLT_TOL     <- 1e-12

  # Main figure text uses a 12-point base size; dense text inside the dose boxes
  # is slightly smaller.
  BASE_PT <- 12

  COL_CONSISTENT    <- "#B65C57"
  COL_INCONSISTENT  <- "#4F78A6"
  FILL_CONSISTENT   <- "#F5E3E1"
  FILL_INCONSISTENT <- "#E3EDF6"
  FILL_NEUTRAL      <- "#FFFFFF"
  COL_BORDER        <- "#555555"
  COL_SEPARATOR     <- "#B8B8B8"
  COL_TEXT          <- "#222222"

  scenario_order <- c(
    "E0", "E1",
    "E2", "E3",
    "M1", "M2",
    "L1", "L2",
    "S1", "S2"
  )

  # Fixed horizontal positions.
  # The larger spacing leaves enough room for "(AUC, DLT probability)" values.
  X_DOSE <- seq(from = 1.25,to = 5.15,length = 5)

  Y_PED   <- 0.82
  Y_ADULT <- 0.25
  Y_DIV   <- 0.50

  # -----------------------------------------------------------------------------
  # 2. Scenario truths
  # -----------------------------------------------------------------------------

  pediatric_auc_base <- c(0.50, 0.70, 1.00, 1.35, 1.80)
  pediatric_dlt_base <- c(0.09, 0.17, 0.30, 0.45, 0.60)
  adult_auc_e1       <- c(0.35, 0.50, 0.70, 1.00, 1.35)

  scenarios <- list(
    E0  = list(
      purpose  = "Exact exposure and toxicity alignment",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = pediatric_dlt_base,
      p_auc    = pediatric_auc_base,
      a_auc    = c(0.50, 0.70, 1.00, 1.35, 1.80)
    ),
    E1  = list(
      purpose  = "Systematic exposure shift",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = c(0.04, 0.09, 0.17, 0.30, 0.45),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    E2  = list(
      purpose  = "Larger exposure shift",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = c(0.02, 0.04, 0.09, 0.17, 0.30),
      p_auc    = pediatric_auc_base,
      a_auc    = c(0.25, 0.35, 0.50, 0.70, 1.00)
    ),
    E3  = list(
      purpose  = "Irregular and many-to-one exposure mapping",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = c(0.09, 0.13, 0.17, 0.37, 0.60),
      p_auc    = pediatric_auc_base,
      a_auc    = c(0.50, 0.60, 0.70, 1.17, 1.80)
    ),
    M1  = list(
      purpose  = "Pediatric MTD at the lowest dose",
      true_mtd = 1L,
      p_dlt    = c(0.30, 0.45, 0.60, 0.72, 0.82),
      a_dlt    = c(0.17, 0.30, 0.45, 0.60, 0.72),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    M2  = list(
      purpose  = "Pediatric MTD at the highest dose",
      true_mtd = 5L,
      p_dlt    = c(0.03, 0.06, 0.10, 0.17, 0.30),
      a_dlt    = c(0.01, 0.03, 0.06, 0.10, 0.17),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    L1  = list(
      purpose  = "Local incompatibility below the target dose",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = c(0.04, 0.09, 0.27, 0.30, 0.45),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    L2  = list(
      purpose  = "Local incompatibility at the target dose",
      true_mtd = 3L,
      p_dlt    = pediatric_dlt_base,
      a_dlt    = c(0.04, 0.09, 0.17, 0.20, 0.45),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    S1  = list(
      purpose  = "No acceptable pediatric dose;\ncompatible matched adult information",
      true_mtd = NA_integer_,
      p_dlt    = c(0.35, 0.45, 0.55, 0.65, 0.75),
      a_dlt    = c(0.25, 0.35, 0.45, 0.55, 0.65),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    ),
    S2  = list(
      purpose  = "No acceptable pediatric dose;\nlower adult toxicity at matched exposures",
      true_mtd = NA_integer_,
      p_dlt    = c(0.35, 0.45, 0.55, 0.65, 0.75),
      a_dlt    = c(0.15, 0.25, 0.35, 0.45, 0.55),
      p_auc    = pediatric_auc_base,
      a_auc    = adult_auc_e1
    )
  )

  # -----------------------------------------------------------------------------
  # 3. Truth-level AUC matching
  # -----------------------------------------------------------------------------

  oracle_auc_match <- function(
    p_auc, adult_auc,
    caliper = AUC_CALIPER,
    tie_tol = AUC_TIE_TOL
  ) {

    dist <- abs(log(p_auc / adult_auc))
    dmin <- min(dist)

    if (dmin > caliper) {
      return(NA_integer_)
    }

    candidates <- which(abs(dist - dmin) <= tie_tol)

    # Tie rule: higher adult AUC, then higher adult nominal dose.
    # Because adult AUC is ordered by dose in these scenarios,
    # the larger index implements both criteria.
    max(candidates)
  }

  is_dlt_consistent <- function(p_dlt, a_dlt, tol  = DLT_TOL) {
    abs(p_dlt - a_dlt) <= tol
  }

  # -----------------------------------------------------------------------------
  # 4. Build plotting data
  # -----------------------------------------------------------------------------

  dose_rows <- list()
  link_rows <- list()
  meta_rows <- list()

  for (sid in scenario_order) {

    s <- scenarios[[sid]]

    for (j in 1:5) {

      dose_rows[[length(dose_rows) + 1L]] <- data.frame(
        scenario         = sid,
        population       = "Pediatric",
        dose_index       = j,
        x                = X_DOSE[j],
        y                = Y_PED,
        auc              = s$p_auc[j],
        dlt              = s$p_dlt[j],
        is_mtd           = !is.na(s$true_mtd) && j == s$true_mtd,
        stringsAsFactors = FALSE
      )

      dose_rows[[length(dose_rows) + 1L]] <- data.frame(
        scenario         = sid,
        population       = "Adult",
        dose_index       = j,
        x                = X_DOSE[j],
        y                = Y_ADULT,
        auc              = s$a_auc[j],
        dlt              = s$a_dlt[j],
        is_mtd           = FALSE,
        stringsAsFactors = FALSE
      )
    }

    matches <- vapply(
      s$p_auc,
      oracle_auc_match,
      integer(1),
      adult_auc  = s$a_auc
    )

    for (j in 1:5) {

      k <- matches[j]

      # No acceptable match -> no arrow.
      if (is.na(k)) next

      status <- if (
        is_dlt_consistent(s$p_dlt[j], s$a_dlt[k])
      ) {
        "AUC matched; DLT consistent"
      } else {
        "AUC matched; DLT inconsistent"
      }

      link_rows[[length(link_rows) + 1L]] <- data.frame(
        scenario         = sid,
        p_index          = j,
        a_index          = k,
        x                = X_DOSE[j],
        xend             = X_DOSE[k],
        y                = 0.64,
        yend             = 0.42,
        status           = status,
        stringsAsFactors = FALSE
      )
    }

    meta_rows[[length(meta_rows) + 1L]] <- data.frame(
      scenario         = sid,
      purpose          = s$purpose,
      stringsAsFactors = FALSE
    )
  }

  dose_df <- do.call(rbind, dose_rows)
  link_df <- do.call(rbind, link_rows)
  meta_df <- do.call(rbind, meta_rows)

  dose_df$scenario <- factor(dose_df$scenario, levels  = scenario_order)
  link_df$scenario <- factor(link_df$scenario, levels  = scenario_order)
  meta_df$scenario <- factor(meta_df$scenario, levels  = scenario_order)

  # -----------------------------------------------------------------------------
  # 5. Adult box colors
  # -----------------------------------------------------------------------------

  dose_df$box_status <- "Neutral"

  for (sid in scenario_order) {

    links_s <- link_df[link_df$scenario == sid, , drop  = FALSE]
    if (nrow(links_s) == 0L) next

    for (k in unique(links_s$a_index)) {

      incoming_status <- unique(
        as.character(links_s$status[links_s$a_index == k])
      )

      # If multiple pediatric doses map to the same adult dose and at least one
      # matched pair is DLT-inconsistent, show the adult box in blue.
      box_status <- if (
        "AUC matched; DLT inconsistent" %in% incoming_status
      ) {
        "AUC matched; DLT inconsistent"
      } else {
        "AUC matched; DLT consistent"
      }

      idx <- dose_df$scenario == sid &
      dose_df$population == "Adult" &
      dose_df$dose_index == k

      dose_df$box_status[idx] <- box_status
    }
  }

  dose_df$box_status <- factor(
    dose_df$box_status,
    levels  = c(
      "Neutral",
      "AUC matched; DLT consistent",
      "AUC matched; DLT inconsistent"
    )
  )

  link_df$status <- factor(
    link_df$status,
    levels  = c(
      "AUC matched; DLT consistent",
      "AUC matched; DLT inconsistent"
    )
  )

  # -----------------------------------------------------------------------------
  # 6. Box labels and geometry
  # -----------------------------------------------------------------------------

  dose_df$dose_name <- ifelse(
    dose_df$population == "Pediatric",
    paste0("P", dose_df$dose_index, ifelse(dose_df$is_mtd, "*", "")),
    paste0("A", dose_df$dose_index)
  )

  # The row label and figure note define the ordered pair shown in each box.
  dose_df$value_label <- sprintf("(%.2f, %.2f)", dose_df$auc, dose_df$dlt)

  BOX_HALF_WIDTH  <- 0.46
  BOX_HALF_HEIGHT <- 0.175

  dose_df$xmin <- dose_df$x - BOX_HALF_WIDTH
  dose_df$xmax <- dose_df$x + BOX_HALF_WIDTH
  dose_df$ymin <- dose_df$y - BOX_HALF_HEIGHT
  dose_df$ymax <- dose_df$y + BOX_HALF_HEIGHT

  # -----------------------------------------------------------------------------
  # 7. Plot
  # -----------------------------------------------------------------------------

  p <- ggplot() +

  geom_hline(
    yintercept = Y_DIV,
    linewidth  = 0.32,
    linetype   = "dashed",
    colour     = COL_SEPARATOR
  ) +

  geom_rect(
    data  = dose_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = box_status
    ),
    colour    = COL_BORDER,
    linewidth = 0.40
  ) +
  
  geom_segment(
    data  = link_df,
    aes(
      x      = x,
      y      = y,
      xend   = xend,
      yend   = yend,
      colour = status
    ),
    linewidth = 0.55,
    lineend   = "round",
    arrow     = grid::arrow(
      length = grid::unit(0.075, "inches"),
      type   = "open"
    )
  ) +
  
  # Dose level
  geom_text(
    data  = dose_df,
    aes(
      x        = x,
      y        = y + 0.066,
      label    = dose_name,
      fontface = ifelse(is_mtd, "bold", "plain")
    ),
    size   = 3.80,
    colour = COL_TEXT
  ) +

  # Ordered pair: (AUC, DLT probability)
  geom_text(
    data  = dose_df,
    aes(
      x     = x,
      y     = y - 0.064,
      label = value_label
    ),
    size   = 3.30,
    colour = COL_TEXT
  ) +

  # Row labels: terminology exactly matches the figure note.
  geom_text(
    data  = meta_df,
    aes(x  = -0.2, y  = Y_PED),
    label      = "Pediatric\n(AUC, DLT)",
    hjust      = 0,
    vjust      = 0.5,
    size       = 3.2,
    fontface   = "bold",
    colour     = COL_TEXT,
    lineheight = 0.95
  ) +

  geom_text(
    data  = meta_df,
    aes(x  = -0.2, y  = Y_ADULT),
    label      = "Adult\n(AUC, DLT)",
    hjust      = 0,
    vjust      = 0.5,
    size       = 3.2,
    fontface   = "bold",
    colour     = COL_TEXT,
    lineheight = 0.95
  ) +

  # Scenario ID
  geom_text(
    data  = meta_df,
    aes(x  = 0.05, y  = 1.22, label  = scenario),
    hjust    = 0,
    size     = 4.55,
    fontface = "bold",
    colour   = COL_TEXT
  ) +

  # Scenario purpose
  geom_text(
    data  = meta_df,
    aes(x  = 3.15, y  = 1.22, label  = purpose),
    hjust      = 0.5,
    size       = 4,
    colour     = COL_TEXT,
    lineheight = 0.95
  ) +

  facet_wrap(
    ~ scenario,
    ncol   = 2,
    scales = "fixed"
  ) +

  scale_fill_manual(
    values  = c(
      "Neutral"  = FILL_NEUTRAL,
      "AUC matched; DLT consistent"  = FILL_CONSISTENT,
      "AUC matched; DLT inconsistent"  = FILL_INCONSISTENT
    ),
    breaks  = c(
      "AUC matched; DLT consistent",
      "AUC matched; DLT inconsistent"
    ),
    labels  = c(
      "AUC matched and DLT probabilities equal",
      "AUC matched but DLT probabilities unequal"
    ),
    name  = NULL
  ) +

  scale_colour_manual(
    values  = c(
      "AUC matched; DLT consistent"  = COL_CONSISTENT,
      "AUC matched; DLT inconsistent"  = COL_INCONSISTENT
    ),
    guide  = "none"
  ) +

  coord_cartesian(
    xlim = c(0.00, 5.38),
    ylim = c(0.00, 1.34),
    clip = "off"
  ) +

  labs(
    caption  = paste0(
      "Note: * True pediatric MTD. \n",
      "           Box values are shown as (AUC, DLT probability). \n",
      "           No connecting arrow indicates no acceptable adult AUC match ",
      "within c = ", sprintf("%.2f", AUC_CALIPER), "."
    )
  ) +

  guides(
    fill  = guide_legend(
      nrow         = 1,
      byrow        = TRUE,
      override.aes = list(
        colour    = COL_BORDER,
        linewidth = 0.4
      )
    )
  ) +

  theme_void(base_size  = BASE_PT) +

  theme(
    text  = element_text(
      family = "sans",
      colour = COL_TEXT,
      size   = BASE_PT
    ),

    strip.text       = element_blank(),
    strip.background = element_blank(),

    panel.border  = element_rect(
      colour    = "#4D4D4D",
      fill      = NA,
      linewidth = 0.52
    ),

    panel.spacing.x = grid::unit(0.055, "inches"),
    panel.spacing.y = grid::unit(0.040, "inches"),

    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.text       = element_text(size  = 11.3),
    legend.key.width  = grid::unit(0.32, "inches"),
    legend.key.height = grid::unit(0.18, "inches"),
    legend.spacing.x  = grid::unit(0.12, "inches"),

    plot.caption  = element_text(
      size       = 10.7,
      hjust      = 0,
      lineheight = 1.05,
      margin     = margin(t  = 7)
    ),

    plot.margin  = margin(
      t    = 7,
      r    = 7,
      b    = 7,
      l    = 7,
      unit = "pt"
    )
  )

  print(p)

  ggplot2::ggsave(
    filename    = output_tiff,
    plot        = p,
    width       = 8.55,
    height      = 9.70,
    units       = "in",
    dpi         = 800,
    device      = "tiff",
    compression = "lzw",
    bg          = "white"
  )
  # 9. QC table
  # -----------------------------------------------------------------------------

  mapping_qc <- do.call(
    rbind,
    lapply(scenario_order, function(sid) {

        s <- scenarios[[sid]]

        k <- vapply(
          s$p_auc,
          oracle_auc_match,
          integer(1),
          adult_auc  = s$a_auc
        )

        data.frame(
          Scenario          = sid,
          Pediatric         = paste0("P", 1:5),
          Adult_match       = ifelse(is.na(k), "None", paste0("A", k)),
          Pediatric_AUC     = s$p_auc,
          Matched_adult_AUC = ifelse(is.na(k), NA_real_, s$a_auc[k]),
          Pediatric_DLT     = s$p_dlt,
          Matched_adult_DLT = ifelse(is.na(k), NA_real_, s$a_dlt[k]),
          DLT_status        = ifelse(
            is.na(k),
            "No match",
            ifelse(
              abs(s$p_dlt - s$a_dlt[k]) <= DLT_TOL,
              "Consistent",
              "Inconsistent"
            )
          ),
          stringsAsFactors  = FALSE
        )
    })
  )

  print(mapping_qc, row.names  = FALSE)

  message("Figure 1 saved to: ", normalizePath(output_tiff, winslash  = "/", mustWork  = FALSE))
  invisible(list(path  = output_tiff, plot  = p, mapping_qc  = mapping_qc))
}

# =============================================================================
# Figure 2: dose-finding performance
# =============================================================================
create_figure2 <- function(analysis, output_tiff) {
  required_pkgs <- c("dplyr", "tidyr", "ggplot2")
  missing_pkgs  <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly  = TRUE)
  ]
  if (length(missing_pkgs) > 0L) {
    stop(
      "Required packages are unavailable: ",
      paste(missing_pkgs, collapse  = ", "),
      call.  = FALSE
    )
  }
  suppressPackageStartupMessages({
      library(dplyr)
      library(tidyr)
      library(ggplot2)
  })
  dir.create(dirname(output_tiff), recursive  = TRUE, showWarnings  = FALSE)
  FIG_WIDTH_MM  <- 180
  FIG_HEIGHT_MM <- 115
  FIG_DPI       <- 800
  BASE_FONT_PT  <- 12

  # Keep the figure at the intended final publication width.
  # If the journal later rescales the artwork, all text and line widths
  # will be rescaled proportionally.


  if (is.null(analysis$main$main)) {
    stop("analysis must be returned by run_formal_full_analysis().", call.  = FALSE)
  }
  main_res       <- analysis$main$main
  scenario_order <- c(
    "E0", "E1", "E2", "E3",
    "M1", "M2",
    "L1", "L2"
  )

  method_order <- c(
    "Keyboard",
    "PA-Keyboard",
    "PA-CRM"
  )

  required_columns <- c(
    "Scenario",
    "Method",
    "PCS (%)",
    "Patients at overly toxic doses (%)"
  )

  missing_columns <- setdiff(required_columns, names(main_res))

  if (length(missing_columns) > 0) {
    stop(
      "The following required columns are missing from 01_Main: ",
      paste(missing_columns, collapse  = ", ")
    )
  }

  # ---- 4. Prepare plotting data ------------------------------------------------

  plot_dat <- main_res %>%
  filter(
    Scenario %in% scenario_order,
    Method %in% method_order
  ) %>%
  transmute(
    Scenario  = factor(
      Scenario,
      levels  = scenario_order
    ),
    Method  = factor(
      Method,
      levels  = method_order
    ),
    PCS  = `PCS (%)`,
    `Overtoxic allocation`  = 
    `Patients at overly toxic doses (%)`
  ) %>%
  pivot_longer(
    cols  = c(
      PCS,
      `Overtoxic allocation`
    ),
    names_to  = "Outcome",
    values_to = "Value"
  ) %>%
  mutate(
    Outcome  = factor(
      Outcome,
      levels  = c(
        "PCS",
        "Overtoxic allocation"
      )
    )
  )

  # ---- 5. Visual encoding ------------------------------------------------------

  # Method is encoded primarily by color and secondarily by point shape.
  # Outcome is encoded by line type; the overtoxic-allocation curves are
  # additionally lightened to reduce visual competition with PCS.

  method_colors <- c(
    "Keyboard"     = "#6B6B6B",
    "PA-Keyboard"  = "#3F6F8C",
    "PA-CRM"       = "#A46A55"
  )

  method_shapes <- c(
    "Keyboard"     = 16,
    "PA-Keyboard"  = 15,
    "PA-CRM"       = 17
  )

  # ---- 6. Main figure ----------------------------------------------------------

  p <- ggplot(
    plot_dat,
    aes(
      x        = Scenario,
      y        = Value,
      group    = interaction(Method, Outcome),
      color    = Method,
      shape    = Method,
      linetype = Outcome,
      alpha    = Outcome
    )
  ) +

  # Light separators between E, M, and L scenario groups
  geom_vline(
    xintercept = c(4.5, 6.5),
    linewidth  = 0.25,
    color      = "grey82"
  ) +

  geom_line(
    linewidth  = 0.60
  ) +

  geom_point(
    size   = 2.45,
    stroke = 0.30
  ) +

  # Group headings. geom_text/annotate uses mm rather than pt,
  # so divide by ggplot2::.pt to obtain the intended point size.
  annotate(
    "text",
    x     = 2.5,
    y     = 92.5,
    label = "Exposure correspondence",
    size  = 11.5 / ggplot2::.pt,
    color = "grey35"
  ) +

  annotate(
    "text",
    x     = 5.5,
    y     = 92.5,
    label = "Boundary MTD",
    size  = 11.5 / ggplot2::.pt,
    color = "grey35"
  ) +

  annotate(
    "text",
    x     = 7.5,
    y     = 92.5,
    label = "Local incompatibility",
    size  = 11.5 / ggplot2::.pt,
    color = "grey35"
  ) +

  scale_color_manual(
    name   = "Method",
    values = method_colors,
    breaks = method_order
  ) +

  scale_shape_manual(
    name   = "Method",
    values = method_shapes,
    breaks = method_order
  ) +

  scale_linetype_manual(
    name   = "Outcome",
    values = c(
      "PCS"  = "solid",
      "Overtoxic allocation"  = "dashed"
    )
  ) +

  scale_alpha_manual(
    name   = "Outcome",
    values = c(
      "PCS"  = 1.00,
      "Overtoxic allocation"  = 0.6
    )
  ) +

  scale_y_continuous(
    name   = "Percentage (%)",
    breaks = seq(0, 90, by  = 10),
    limits = c(0, 96),
    expand = expansion(mult  = c(0, 0))
  ) +

  labs(
    x  = "Scenario"
  ) +

  # Two concise legend rows:
  # row 1 = method; row 2 = outcome
  guides(
    color  = guide_legend(
      order        = 1,
      nrow         = 1,
      override.aes = list(
        linetype  = "solid",
        alpha     = 1,
        linewidth = 0.60
      )
    ),

    shape  = guide_legend(
      order = 1,
      nrow  = 1
    ),

    linetype  = guide_legend(
      order        = 2,
      nrow         = 1,
      override.aes = list(
        color     = "grey30",
        shape     = NA,
        linewidth = 0.60
      )
    ),

    alpha  = guide_legend(
      order = 2,
      nrow  = 1
    )
  ) +

  # 12 pt is the target base size at the FINAL 180-mm canvas.
  theme_classic(
    base_size   = BASE_FONT_PT,
    base_family = "sans"
  ) +

  theme(
    # Axis titles and tick labels
    axis.title  = element_text(
      size  = BASE_FONT_PT,
      color = "black"
    ),

    axis.text  = element_text(
      size  = BASE_FONT_PT,
      color = "black"
    ),

    axis.ticks  = element_line(
      linewidth = 0.35,
      color     = "black"
    ),

    axis.ticks.length  = unit(
      1.7,
      "mm"
    ),

    axis.line  = element_line(
      linewidth = 0.35,
      color     = "black"
    ),

    # Legends
    legend.position = "bottom",
    legend.box      = "vertical",

    legend.title  = element_text(
      size  = BASE_FONT_PT
    ),

    legend.text  = element_text(
      size  = BASE_FONT_PT
    ),

    legend.key.width  = unit(
      10,
      "mm"
    ),

    legend.key.height  = unit(
      4.5,
      "mm"
    ),

    legend.spacing.y  = unit(
      1.0,
      "mm"
    ),

    legend.margin  = margin(
      t    = 1.5,
      r    = 0,
      b    = 0,
      l    = 0,
      unit = "mm"
    ),

    # Compact journal-style margins
    plot.margin  = margin(
      t    = 3.0,
      r    = 3.0,
      b    = 2.0,
      l    = 2.5,
      unit = "mm"
    )
  )

  # ---- 7. Display --------------------------------------------------------------

  print(p)
  ggsave(
    filename    = output_tiff,
    plot        = p,
    width       = FIG_WIDTH_MM,
    height      = FIG_HEIGHT_MM,
    units       = "mm",
    dpi         = FIG_DPI,
    device      = "tiff",
    compression = "lzw",
    bg          = "white"
  )

  message("Figure 2 saved to: ", normalizePath(output_tiff, winslash  = "/", mustWork  = FALSE))
  invisible(list(path  = output_tiff, plot  = p, data  = plot_dat))
}

