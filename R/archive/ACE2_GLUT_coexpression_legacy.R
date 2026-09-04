# ─────────────────────────────────────────────────────────────────────────────
# ACE2 – GLUT co-expression scatter plots
# Adjust ALL visual options in the CONFIG block below, then run the full script
# ─────────────────────────────────────────────────────────────────────────────

library(readxl); library(ggplot2); library(janitor); library(dplyr)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIG — change anything here, leave the plotting code untouched
# ══════════════════════════════════════════════════════════════════════════════

DATA_FILE   <- "GDM_aspirin_data.xlsx"   # path to your Excel file
SHEET       <- "Sheet1"                  # sheet name
SKIP_ROWS   <- 1                         # rows to skip before header

# ── Output ────────────────────────────────────────────────────────────────────
OUT_A       <- "Figure4A_ACE2_GLUT1.png"
OUT_B       <- "Figure4B_ACE2_GLUT3.png"
FIG_WIDTH   <- 6     # inches
FIG_HEIGHT  <- 5.0     # inches
DPI         <- 300

# ── Text sizes (points) ───────────────────────────────────────────────────────
SIZE_BASE        <- 12   # base size — all other sizes scale from this
SIZE_TITLE       <- 12   # plot title
SIZE_AXIS_TITLE  <- 11   # x and y axis labels
SIZE_AXIS_TEXT   <- 10   # tick mark numbers
SIZE_LEGEND_TITLE<- 10   # "Group" / "Aspirin" legend headers
SIZE_LEGEND_TEXT <-  9   # legend item labels
SIZE_ANNOT       <-  4   # r and p annotation inside plot (in ggplot units)

# ── Colors (Okabe-Ito colorblind-safe palette) ────────────────────────────────
COL_CONTROL  <- "#56B4E9"   # sky blue
COL_DIET     <- "#E69F00"   # amber
COL_INSULIN  <- "#CC79A7"   # mauve/pink

# ── Point appearance ──────────────────────────────────────────────────────────
SHAPE_NO_ASP <- 16     # circle  = no aspirin
SHAPE_ASP    <- 17     # triangle = aspirin exposed
POINT_SIZE   <- 3.5    # point diameter
POINT_ALPHA  <- 0.92   # transparency (0 = invisible, 1 = solid)
POINT_STROKE <- 0.7    # white edge thickness around each point

# ── Regression line ───────────────────────────────────────────────────────────
LINE_COLOR   <- "grey35"
LINE_TYPE    <- "dashed"   # "solid" | "dashed" | "dotted"
LINE_WIDTH   <- 0.9
CI_FILL      <- "grey80"   # 95% confidence band fill
CI_ALPHA     <- 0.25

# ── Axis labels ───────────────────────────────────────────────────────────────
XLAB <- "ACE2 protein (ELISA, pg/mg total protein)"
YLAB_A <- "GLUT1 protein (ELISA, pg/mg total protein)"
YLAB_B <- "GLUT3 protein (ELISA, pg/mg total protein)"

# ── Plot titles ───────────────────────────────────────────────────────────────
TITLE_A <- "A   ACE2 – GLUT1 co-expression"
TITLE_B <- "B   ACE2 – GLUT3 co-expression"

# ── Legend position ───────────────────────────────────────────────────────────
LEGEND_POS <- "right"   # "right" | "bottom" | "none"

# ══════════════════════════════════════════════════════════════════════════════
# DATA LOADING  
# ══════════════════════════════════════════════════════════════════════════════

raw <- read_excel(DATA_FILE, sheet = SHEET, skip = SKIP_ROWS)
colnames(raw) <- c("ID","Notes","Group","Aspirin","GestAge",
                   "PrePregBMI","CurPregBMI",
                   "ACE2_ELISA","GLUT1_ELISA","GLUT3_ELISA",
                   "GLUT4_ELISA","GLUT9_ELISA",
                   "ACE2_WB","GLUT1_WB","GLUT3_WB","GLUT4_WB","GLUT9_WB")

df <- raw |>
  
  # Keep only rows where ID is a plain integer — drops header and blank rows
  filter(grepl("^\\d+$", as.character(ID))) |>
  mutate(
    across(c(ACE2_ELISA, GLUT1_ELISA, GLUT3_ELISA,
             GLUT4_ELISA, GLUT9_ELISA,
             ACE2_WB, GLUT1_WB, GLUT3_WB, GLUT4_WB, GLUT9_WB),
           as.numeric),
    Group   = factor(Group,   levels = c("Control","GDM_diet","GDM_insulin")),
    Aspirin = factor(Aspirin, levels = c("No","Yes"))
  ) |>
  filter(!is.na(ACE2_ELISA), !is.na(GLUT1_ELISA), !is.na(GLUT3_ELISA),
         !is.na(Group), !is.na(Aspirin))
# ══════════════════════════════════════════════════════════════════════════════
# PLOT FUNCTION
# ══════════════════════════════════════════════════════════════════════════════

make_scatter <- function(data, y_var, y_label, plot_title) {
  
  # Pearson r and p
  ct    <- cor.test(data$ACE2_ELISA, data[[y_var]], method = "pearson")
  r_val <- round(ct$estimate, 2)
  p_val <- ct$p.value
  p_str <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  ann   <- sprintf("italic(r) == %s * ', ' * '%s'", r_val, p_str)
  
  ggplot(data, aes(x = ACE2_ELISA, y = .data[[y_var]])) +
    
    # Regression line + CI band
    geom_smooth(
      method    = "lm",
      se        = TRUE,
      color     = LINE_COLOR,
      fill      = CI_FILL,
      alpha     = CI_ALPHA,
      linetype  = LINE_TYPE,
      linewidth = LINE_WIDTH,
      show.legend = FALSE
    ) +
    
    # Data points
    geom_point(
      aes(color = Group, shape = Aspirin),
      size   = POINT_SIZE,
      stroke = POINT_STROKE,
      alpha  = POINT_ALPHA
    ) +
    
    # r and p annotation
    annotate(
      "text",
      x     = Inf, y = -Inf,
      hjust = 1.08, vjust = -0.7,
      label = ann,
      parse = TRUE,
      size  = SIZE_ANNOT,
      fontface = "italic"
    ) +
    
    # Color scale
    scale_color_manual(
      na.translate = FALSE,          
      values = c(
        Control     = COL_CONTROL,
        GDM_diet    = COL_DIET,
        GDM_insulin = COL_INSULIN
      ),
      labels = c(
        Control     = "Control (n = 3)",
        GDM_diet    = "GDM \u2013 diet (n = 7)",
        GDM_insulin = "GDM \u2013 insulin (n = 8)"
      )
    ) +
    
    # Shape scale
    scale_shape_manual(
      na.translate = FALSE,          
      values = c(No = SHAPE_NO_ASP, Yes = SHAPE_ASP),
      labels = c(No = "No aspirin", Yes = "Aspirin exposed")
    ) +
    
    # Labels
    labs(
      title  = plot_title,
      x      = XLAB,
      y      = y_label,
      color  = "Group",
      shape  = "Aspirin"
    ) +
    
    # Theme
    theme_classic(base_size = SIZE_BASE) +
    theme(
      plot.title     = element_text(face = "bold",  size = SIZE_TITLE,        hjust = 0),
      axis.title     = element_text(size = SIZE_AXIS_TITLE),
      axis.text      = element_text(size = SIZE_AXIS_TEXT),
      legend.title   = element_text(size = SIZE_LEGEND_TITLE, face = "bold"),
      legend.text    = element_text(size = SIZE_LEGEND_TEXT),
      legend.position = LEGEND_POS,
      axis.line      = element_line(linewidth = 0.4),
      axis.ticks     = element_line(linewidth = 0.4)
    )
}

# ══════════════════════════════════════════════════════════════════════════════
# SAVE PLOTS
# ══════════════════════════════════════════════════════════════════════════════

pA <- make_scatter(df, "GLUT1_ELISA", YLAB_A, TITLE_A)
pB <- make_scatter(df, "GLUT3_ELISA", YLAB_B, TITLE_B)

ggsave(OUT_A, plot = pA, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = DPI)
ggsave(OUT_B, plot = pB, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = DPI)

message("Saved: ", OUT_A)
message("Saved: ", OUT_B)
