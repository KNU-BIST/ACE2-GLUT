# =============================================================================
# Phase 7 — Clinical placental ACE2–GLUT protein correlations
# Script: 07_clinical_ACE2_GLUT_correlations.R
# =============================================================================
#
# Project
# -------
# Placental ACE2 and glucose transporter expression in gestational diabetes
# mellitus
#
# Author
# ------
# Md. Masudul Haque
#
# Purpose
# -------
# Generate patient-level ACE2–GLUT1 and ACE2–GLUT3 protein correlation plots
# from the institutional placental cohort, calculate full-cohort Pearson
# correlations, and reproduce the sensitivity/within-group analyses used to
# evaluate between-group clustering.
#
# Reproducibility boundary
# ------------------------
# The biological unit of analysis is one placenta per pregnancy.
#
# Technical replicates from ELISA or other laboratory measurements must be
# summarized at the placenta level before this script is run and are not
# treated as independent biological observations.
#
# These analyses are observational and exploratory. Correlation does not
# establish regulatory direction or a causal ACE2–GLUT relationship.
#
# Aspirin exposure is represented visually but is not interpreted as a
# randomized treatment effect.
#
# Inputs
# ------
# data/private/clinical/GDM_aspirin_protein_data.xlsx
#
# Outputs
# -------
# figures/main/Figure4A_ACE2_GLUT1.png
# figures/main/Figure4B_ACE2_GLUT3.png
# outputs/statistics/clinical/ACE2_GLUT_full_correlations.csv
# outputs/statistics/clinical/ACE2_GLUT_sensitivity_within_group.csv
#
# Cohort
# ------
# Total n = 18 independent placentas
#   Control      = 3
#   GDM-diet     = 7
#   GDM-insulin  = 8
#
# =============================================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Install package 'here' before running this script.")
}
source(here::here("R", "00_project_setup.R"))

data_file <- file.path(
  path_data_private,
  "GDM_aspirin_protein_data.xlsx"
)

assert_file_exists(data_file)

sheet_name <- "Sheet1"
skip_rows <- 1

# -----------------------------------------------------------------------------
# Figure configuration
# -----------------------------------------------------------------------------

fig_width <- 6
fig_height <- 5

size_base <- 12
size_title <- 12
size_axis_title <- 11
size_axis_text <- 10
size_legend_title <- 10
size_legend_text <- 9
size_annot <- 4

col_control <- "#56B4E9"
col_diet <- "#E69F00"
col_insulin <- "#CC79A7"

shape_no_asp <- 16
shape_asp <- 17

point_size_clinical <- 3.5
point_alpha <- 0.92
point_stroke <- 0.7

line_color <- "grey35"
line_type <- "dashed"
line_width <- 0.9
ci_fill <- "grey80"
ci_alpha <- 0.25

x_label <- "ACE2 protein (ELISA, pg/mg total protein)"
y_label_glut1 <- "GLUT1 protein (ELISA, pg/mg total protein)"
y_label_glut3 <- "GLUT3 protein (ELISA, pg/mg total protein)"

# -----------------------------------------------------------------------------
# 1. Load and validate the institutional protein dataset
# -----------------------------------------------------------------------------

raw <- readxl::read_excel(
  data_file,
  sheet = sheet_name,
  skip = skip_rows
)

expected_columns <- c(
  "ID",
  "Notes",
  "Group",
  "Aspirin",
  "GestAge",
  "PrePregBMI",
  "CurPregBMI",
  "ACE2_ELISA",
  "GLUT1_ELISA",
  "GLUT3_ELISA",
  "GLUT4_ELISA",
  "GLUT9_ELISA",
  "ACE2_WB",
  "GLUT1_WB",
  "GLUT3_WB",
  "GLUT4_WB",
  "GLUT9_WB"
)

if (ncol(raw) != length(expected_columns)) {
  stop(
    "Unexpected number of columns in clinical protein spreadsheet. Expected ",
    length(expected_columns),
    " but found ",
    ncol(raw),
    "."
  )
}

colnames(raw) <- expected_columns

protein_columns <- c(
  "ACE2_ELISA",
  "GLUT1_ELISA",
  "GLUT3_ELISA",
  "GLUT4_ELISA",
  "GLUT9_ELISA",
  "ACE2_WB",
  "GLUT1_WB",
  "GLUT3_WB",
  "GLUT4_WB",
  "GLUT9_WB"
)

df <- raw |>
  dplyr::filter(
    grepl(
      "^\\d+$",
      as.character(ID)
    )
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(protein_columns),
      as.numeric
    ),
    Group = factor(
      Group,
      levels = c(
        "Control",
        "GDM_diet",
        "GDM_insulin"
      )
    ),
    Aspirin = factor(
      Aspirin,
      levels = c(
        "No",
        "Yes"
      )
    )
  ) |>
  dplyr::filter(
    !is.na(ACE2_ELISA),
    !is.na(GLUT1_ELISA),
    !is.na(GLUT3_ELISA),
    !is.na(Group),
    !is.na(Aspirin)
  )

if (nrow(df) != 18) {
  warning(
    "Expected 18 biological placental observations but found ",
    nrow(df),
    ". Verify the input spreadsheet before manuscript use."
  )
}

print(
  table(df$Group, useNA = "ifany")
)

# -----------------------------------------------------------------------------
# 2. Full-cohort correlations
# -----------------------------------------------------------------------------

calculate_correlation <- function(
    data,
    x_var,
    y_var,
    analysis_name
) {
  complete <- data[
    complete.cases(
      data[
        ,
        c(
          x_var,
          y_var
        )
      ]
    ),
    ,
    drop = FALSE
  ]

  if (nrow(complete) < 4) {
    return(
      data.frame(
        Analysis = analysis_name,
        Comparison = paste(
          x_var,
          "vs",
          y_var
        ),
        n = nrow(complete),
        r = NA_real_,
        p = NA_real_,
        note = "Too few observations for a stable reported correlation"
      )
    )
  }

  ct <- stats::cor.test(
    complete[[x_var]],
    complete[[y_var]],
    method = "pearson"
  )

  data.frame(
    Analysis = analysis_name,
    Comparison = paste(
      x_var,
      "vs",
      y_var
    ),
    n = nrow(complete),
    r = unname(ct$estimate),
    p = ct$p.value,
    note = ""
  )
}

full_correlations <- dplyr::bind_rows(
  calculate_correlation(
    df,
    "ACE2_ELISA",
    "GLUT1_ELISA",
    "Full cohort"
  ),
  calculate_correlation(
    df,
    "ACE2_ELISA",
    "GLUT3_ELISA",
    "Full cohort"
  ),
  calculate_correlation(
    df,
    "ACE2_ELISA",
    "GLUT4_ELISA",
    "Full cohort"
  ),
  calculate_correlation(
    df,
    "ACE2_ELISA",
    "GLUT9_ELISA",
    "Full cohort"
  )
)

write.csv(
  full_correlations,
  file.path(
    path_stats_clinical,
    "ACE2_GLUT_full_correlations.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3. Sensitivity and within-group correlations
# -----------------------------------------------------------------------------

sensitivity_data <- df |>
  dplyr::filter(
    Group != "GDM_diet"
  )

sensitivity_results <- dplyr::bind_rows(
  calculate_correlation(
    sensitivity_data,
    "ACE2_ELISA",
    "GLUT1_ELISA",
    "Sensitivity: GDM_diet excluded"
  ),
  calculate_correlation(
    sensitivity_data,
    "ACE2_ELISA",
    "GLUT3_ELISA",
    "Sensitivity: GDM_diet excluded"
  )
)

within_group_results <- dplyr::bind_rows(
  lapply(
    levels(df$Group),
    function(group_name) {
      group_df <- df |>
        dplyr::filter(
          Group == group_name
        )

      dplyr::bind_rows(
        calculate_correlation(
          group_df,
          "ACE2_ELISA",
          "GLUT1_ELISA",
          paste0(
            "Within-group: ",
            group_name
          )
        ),
        calculate_correlation(
          group_df,
          "ACE2_ELISA",
          "GLUT3_ELISA",
          paste0(
            "Within-group: ",
            group_name
          )
        )
      )
    }
  )
)

sensitivity_within <- dplyr::bind_rows(
  sensitivity_results,
  within_group_results
)

write.csv(
  sensitivity_within,
  file.path(
    path_stats_clinical,
    "ACE2_GLUT_sensitivity_within_group.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4. Manuscript Figure 4 scatter plots
# -----------------------------------------------------------------------------

make_scatter <- function(
    data,
    y_var,
    y_label,
    plot_title
) {
  ct <- stats::cor.test(
    data$ACE2_ELISA,
    data[[y_var]],
    method = "pearson"
  )

  r_val <- round(
    unname(ct$estimate),
    2
  )

  p_val <- ct$p.value

  p_str <- if (p_val < 0.001) {
    "p < 0.001"
  } else {
    sprintf(
      "p = %.3f",
      p_val
    )
  }

  annotation <- sprintf(
    "italic(r) == %s * ', ' * '%s'",
    r_val,
    p_str
  )

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = ACE2_ELISA,
      y = .data[[y_var]]
    )
  ) +
    ggplot2::geom_smooth(
      method = "lm",
      se = TRUE,
      color = line_color,
      fill = ci_fill,
      alpha = ci_alpha,
      linetype = line_type,
      linewidth = line_width,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = Group,
        shape = Aspirin
      ),
      size = point_size_clinical,
      stroke = point_stroke,
      alpha = point_alpha
    ) +
    ggplot2::annotate(
      "text",
      x = Inf,
      y = -Inf,
      hjust = 1.08,
      vjust = -0.7,
      label = annotation,
      parse = TRUE,
      size = size_annot
    ) +
    ggplot2::scale_color_manual(
      na.translate = FALSE,
      values = c(
        Control = col_control,
        GDM_diet = col_diet,
        GDM_insulin = col_insulin
      ),
      labels = c(
        Control = "Control (n = 3)",
        GDM_diet = "GDM — diet (n = 7)",
        GDM_insulin = "GDM — insulin (n = 8)"
      )
    ) +
    ggplot2::scale_shape_manual(
      na.translate = FALSE,
      values = c(
        No = shape_no_asp,
        Yes = shape_asp
      ),
      labels = c(
        No = "No aspirin",
        Yes = "Aspirin exposed"
      )
    ) +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      color = "Group",
      shape = "Aspirin"
    ) +
    ggplot2::theme_classic(
      base_size = size_base
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = size_title,
        hjust = 0
      ),
      axis.title = ggplot2::element_text(
        size = size_axis_title
      ),
      axis.text = ggplot2::element_text(
        size = size_axis_text
      ),
      legend.title = ggplot2::element_text(
        size = size_legend_title,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = size_legend_text
      ),
      legend.position = "right",
      axis.line = ggplot2::element_line(
        linewidth = 0.4
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.4
      )
    )
}

pA <- make_scatter(
  df,
  "GLUT1_ELISA",
  y_label_glut1,
  "A   ACE2–GLUT1 co-expression"
)

pB <- make_scatter(
  df,
  "GLUT3_ELISA",
  y_label_glut3,
  "B   ACE2–GLUT3 co-expression"
)

save_plot(
  pA,
  file.path(
    path_fig_main,
    "Figure4A_ACE2_GLUT1.png"
  ),
  width = fig_width,
  height = fig_height
)

save_plot(
  pB,
  file.path(
    path_fig_main,
    "Figure4B_ACE2_GLUT3.png"
  ),
  width = fig_width,
  height = fig_height
)

write_session_info("07_clinical_ACE2_GLUT_correlations.R")

message("\nPhase 7 complete.\n")
