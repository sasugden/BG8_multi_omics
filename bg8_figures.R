# Define plot theme ####
plot_theme <- theme(plot.tag = element_text(color="black", face="bold"),
                    strip.text = element_text(size=9),
                    panel.grid = element_blank(),
                    strip.background=element_rect(color="black", size=1),
                    panel.border = element_rect(color="black", fill=NA),
                    axis.text = element_text(color= "black", size=8),
                    axis.title = element_text(color="black", size=9),
                    legend.title = element_text(color="black", size=9),
                    plot.title=element_text(color="black", size=10),
                    legend.text = element_text(color="black", size=8, hjust=0))

#### FIGURE 1: GROWTH CURVES ####
od_data <- read.csv("raw_data/BG8_growth_data.csv")
od_data$Treatment <- factor(od_data$Treatment, levels=c("CH4-AMS", "CH4-NMS", "MeOH-AMS", "MeOH-NMS"))
od_data$Omics <- factor(od_data$Omics, levels=c("Transcriptome", "Metabolome"))

# Define treatment order for plotting ####
od_data <- subset(od_data, is.na(SD)=="FALSE")

vnames <-list(
  'CH4-AMS' = bquote(paste("C", H[4], ":AMS")),
  'CH4-NMS' = bquote(paste("C", H[4], ":NMS")),
  'MeOH-AMS' = bquote(paste("C", H[3], "OH:AMS")),
  'MeOH-NMS' = bquote(paste("C", H[3], "OH:NMS")))

vlabeller <- function(variable,value){
  return(vnames[value])
}

# Plot - with normal and log-transformed axes ####
od_figure <- ggplot(od_data, aes(x=Timepoint, y=Value, color=Omics)) +
  geom_point(size=1.5) +
  geom_line(aes(x=Timepoint, y=Value, color=Omics)) +
  geom_errorbar(aes(x=Timepoint, ymin=Value-SD, ymax=Value+SD), width=3) +
  scale_color_manual(values=c("firebrick1","dodgerblue4"), name="Experiment") +
  scale_y_log10() +
  facet_wrap(~Treatment, labeller=vlabeller) +
  theme_bw() + plot_theme +
  labs(x="\nTime (h)", y="log(OD540)\n") +
  theme(legend.position="right",
        strip.background=element_rect(fill="beige"),
        strip.text=element_text(color=c("black")),
        panel.grid=element_blank())

ggsave("~/od_panels.jpg", od_figure,
       width=6.5, height=5, units="in", dpi=300)

od_data_2 <- od_data
od_data$SD[od_data$Treatment=="MeOH-NMS" & od_data$Omics=="Transcriptome" & od_data$Timepoint==60] <- 0.110
od_data$SD[od_data$Treatment=="MeOH-AMS" & od_data$Omics=="Metabolome" & od_data$Timepoint==48] <- 0.03388

FIGURE_od_log <- ggplot(od_data_2, aes(x=Timepoint, y=Value, color=Omics)) +
  geom_point(size=1.5) +
  #geom_line(data = od_spline, aes(x=Timepoint, y=Value, color=Omics)) +
  ggalt::geom_xspline() +
  geom_errorbar(aes(x=Timepoint, ymin=Value-SD, ymax=Value+SD), width=3) +
  scale_color_manual(values=c("firebrick1","dodgerblue4"), name="Experiment") +
  scale_y_log10() +
  facet_wrap(~Treatment, labeller=vlabeller) +
  theme_bw() +
  labs(x="\nTime (h)", y="log(OD540)\n") +
  theme(axis.title=element_text(size=10, color="black"),
        axis.text=element_text(size=9, color="black"),
        legend.title=element_text(size=10, color="black", hjust=0.5),
        legend.text=element_text(size=9, color="black"),
        legend.text.align=0,
        legend.position="right",
        plot.tag=element_text(face="bold"),
        plot.title=element_text(size=11),
        strip.background=element_rect(fill="beige"),
        strip.text=element_text(color=c("black"), size=10),
        panel.grid=element_blank())

cowplot::save_plot("~/BG8_paper/FIGURE_od_log.jpg", FIGURE_od_log,
       base_width=6.5, base_height=5, units="in", dpi=300)

#### FIGURE 2: METABOLOME + TRANSCRIPTOME PCA SIDE-BY-SIDE ####
pca1 <- ggplot() +
  geom_point(data=rnaseq_pca_results$scores, aes(x=PC1, y=PC2, color=Treatment), size=2) +
  geom_path(data=rnaseq_pca_results$ellipse, aes(x=PC1, y=PC2, color=Treatment)) +
  labs(title = "Transcriptome", tag="a",
       x="PC1 (36.7%)", y="PC2 (21.4%)") +
  scale_color_manual(breaks=levels(rnaseq_sample_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme_bw() + plot_theme +
  theme(legend.position="none")

pca2 <- ggplot() +
  geom_point(data=metab_pca_results$scores, aes(x=PC1, y=PC2, color=Treatment), size=2) +
  geom_path(data=metab_pca_results$ellipse, aes(x=PC1, y=PC2, color=Treatment)) +
  labs(title = "Metabolome", tag="b",
       x="PC1 (68.8%)", y="PC2 (14.1%)") +
  scale_color_manual(breaks=levels(rnaseq_sample_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme_bw() + plot_theme +
  theme(legend.position="none")

pca_legend <- cowplot::get_legend(ggplot() +
  geom_point(data=metab_pca_results$scores, aes(x=PC1, y=PC2, color=Treatment), size=2) +
  geom_path(data=metab_pca_results$ellipse, aes(x=PC1, y=PC2, color=Treatment)) +
  labs(title = "Metabolome", tag="b",
       x="PC1 (68.8%)", y="PC2 (14.1%)") +
  scale_color_manual(breaks=levels(rnaseq_sample_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme_bw() + plot_theme +
  theme(legend.position="right",
        legend.text=element_text(hjust=0)))

grid.arrange(pca1, pca2, pca_legend, nrow=1, widths=c(0.4,0.4,0.2))

FIGURE_pca <- arrangeGrob(pca1, pca2, pca_legend, nrow=1, widths=c(0.4, 0.4, 0.2))
ggsave("FIGURE_pca.jpg", FIGURE_pca, dpi=300, width=6.5, height=3, units="in")

rm(pca1, pca2, pca_legend, FIGURE_pca)

#### FIGURE 3: METHANE VS METHANOL ####
# (a, d) STACKED VOLCANO ####
# Create column identifying significantly differentially abundant items in red.
rnaseq_plot_data <- rnaseq_diff_abund
for(i in c(1:4)){
  rnaseq_plot_data[[i]]$color <- rep("grey")
  for(p in 1:nrow(rnaseq_plot_data[[i]])){
    if(rnaseq_plot_data[[i]][p,"p.adj"] < 0.01 & abs(rnaseq_plot_data[[i]][p,"logFC"]) > 1){
      rnaseq_plot_data[[i]][p,"color"] <- "red"
    }
  }
}


metab_plot_data <- metab_diff_abund
for(i in c(1:4)){
  metab_plot_data[[i]]$color <- rep("grey")
  for(p in 1:nrow(metab_plot_data[[i]])){
    if(metab_plot_data[[i]][p,"p.adj"] < 0.01 & abs(metab_plot_data[[i]][p,"logFC"]) > 1){
      metab_plot_data[[i]][p,"color"] <- "red"
    }
  }
}

# Transcriptome: Methane vs. methanol, NMS
plot1a <- ggplot() +
  geom_point(data=rnaseq_plot_data$CH4vCH3OH_NMS, aes(x=logFC, y=log10q), color=rnaseq_plot_data[["CH4vCH3OH_NMS"]]$color, size=1) +
  #geom_point(data=rnaseq_plot_data$CH4vCH3OH_NMS, aes(x=logFC, y=-1*log10q), color=rnaseq_plot_data[["CH4vCH3OH_NMS"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(y="-log10 q value", tag="a") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,55), expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9, color="white"))

# Metabolome: Methane vs. methanol, NMS
plot1b <- ggplot() +
  #geom_point(data=rnaseq_diff_abund$CH4vCH3OH_NMS, aes(x=logFC, y=log10q), color=rnaseq_diff_abund[["CH4vCH3OH_NMS"]]$color, size=1) +
  geom_point(data=metab_plot_data$CH4vCH3OH_NMS, aes(x=logFC, y=log10q), color=metab_plot_data[["CH4vCH3OH_NMS"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="a") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", color="white"),
        axis.text=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="white"))

plot1a <- ggplotGrob(plot1a)
plot1b <- ggplotGrob(plot1b)
maxWidth = grid::unit.pmax(plot1a$widths[2:5], plot1b$widths[2:5])
plot1a$widths[2:5] <- as.list(maxWidth)
plot1b$widths[2:5] <- as.list(maxWidth)

plot1 <- arrangeGrob(plot1a, plot1b, nrow=2, heights=c(0.45,0.55))
rm(plot1a, plot1b)

# Transcriptome: Methane vs methanol, in AMS
plot4a <- ggplot() +
  geom_point(data=rnaseq_plot_data$CH4vCH3OH_AMS, aes(x=logFC, y=log10q), color=rnaseq_plot_data[["CH4vCH3OH_AMS"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(y="-log10 q value", tag="d") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,55), expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9, color="white"))

plot4b <- ggplot() +
  #geom_point(data=rnaseq_diff_abund$CH4vCH3OH_NMS, aes(x=logFC, y=log10q), color=rnaseq_diff_abund[["CH4vCH3OH_NMS"]]$color, size=1) +
  geom_point(data=metab_plot_data$CH4vCH3OH_AMS, aes(x=logFC, y=log10q), color=metab_plot_data[["CH4vCH3OH_AMS"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="d") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", color="white"),
        axis.text=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="white"))

plot4a <- ggplotGrob(plot4a)
plot4b <- ggplotGrob(plot4b)
maxWidth = grid::unit.pmax(plot4a$widths[2:5], plot4b$widths[2:5])
plot4a$widths[2:5] <- as.list(maxWidth)
plot4b$widths[2:5] <- as.list(maxWidth)

plot4 <- arrangeGrob(plot4a, plot4b, nrow=2, heights=c(0.45,0.55))
rm(plot4a, plot4b)

# (b) RNA COGS - CH4 VS MeOH - NMS ####
plot_data <- rnaseq_numbers$CH4vCH3OH_NMS
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("CH3OH", "Fill", "CH4"))

plot2 <- ggplot() +
  geom_bar(data=plot_data, aes(x=COG, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$COG))) +
  scale_y_continuous(limits=c(0,25), breaks=c(0,6.5,12.5,18.5,25), labels=c(25,18.5,12.5,6.5,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,6.5,12.5,18.5,25))) +
  theme_bw() + plot_theme +
  labs(y="% of DEGs", x="COG", tag="b") +
  scale_fill_manual(values=c("#2e74b6", "white", "#70AD47"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkgreen"),
        axis.ticks.x.top=element_line(color="darkgreen"),
        axis.line.x.top=element_line(color="darkgreen"),
        axis.text.x.bottom=element_text(color="#2e74b6"),
        axis.ticks.x.bottom=element_line(color="#2e74b6"),
        axis.line.x.bottom=element_line(color="#2e74b6"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (e) RNA COGS - CH4 VS MeOH - AMS ####
plot_data <- rnaseq_numbers$CH4vCH3OH_AMS
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("CH3OH", "Fill", "CH4"))

plot5 <- ggplot() +
  geom_bar(data=plot_data, aes(x=COG, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$COG))) +
  scale_y_continuous(limits=c(0,25), breaks=c(0,6.5,12.5,18.5,25), labels=c(25,18.5,12.5,6.5,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,6.5,12.5,18.5,25))) +
  theme_bw() + plot_theme +
  labs(y="% of DEGs", x="COG", tag="e") +
  scale_fill_manual(values=c("#FFC000", "white", "#ED7D31"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="orangered3"),
        axis.ticks.x.top=element_line(color="orangered3"),
        axis.line.x.top=element_line(color="orangered3"),
        axis.text.x.bottom=element_text(color="darkgoldenrod3"),
        axis.ticks.x.bottom=element_line(color="darkgoldenrod3"),
        axis.line.x.bottom=element_line(color="darkgoldenrod3"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (c) METAB COGS - CH4 VS MeOH - NMS ####
plot_data <- metab_numbers$CH4vCH3OH_NMS
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("CH3OH", "Fill", "CH4"))

plot3 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="c") +
  scale_fill_manual(values=c("#2e74b6", "white", "#70AD47"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkgreen"),
        axis.ticks.x.top=element_line(color="darkgreen"),
        axis.line.x.top=element_line(color="darkgreen"),
        axis.text.x.bottom=element_text(color="#2e74b6"),
        axis.ticks.x.bottom=element_line(color="#2e74b6"),
        axis.line.x.bottom=element_line(color="#2e74b6"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (f) METAB COGS - CH4 VS MeOH - AMS ####
plot_data <- metab_numbers$CH4vCH3OH_AMS
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("CH3OH", "Fill", "CH4"))

plot6 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="f") +
  scale_fill_manual(values=c("#FFC000", "white", "#ED7D31"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="orangered3"),
        axis.ticks.x.top=element_line(color="orangered3"),
        axis.line.x.top=element_line(color="orangered3"),
        axis.text.x.bottom=element_text(color="darkgoldenrod3"),
        axis.ticks.x.bottom=element_line(color="darkgoldenrod3"),
        axis.line.x.bottom=element_line(color="darkgoldenrod3"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# Create a legend ####
legend_data <- rbind(metab_numbers$AMSvNMS_Methane, metab_numbers$CH4vCH3OH_AMS)
legend_data <- subset(legend_data, Treatment !="Fill")
legend_data$Treatment <- factor(legend_data$Treatment)

ggplot(legend_data, aes(x=Superfamily, y=Percentage, fill=Treatment)) + geom_bar(stat="identity") +
  scale_fill_manual(breaks=levels(legend_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6"),
                    guide=guide_legend(ncol=2, title.position="top")) +
  theme(legend.position="bottom",
        legend.title=element_text(size=10, color="black", hjust=0.5),
        legend.text=element_text(size=9, color="black"),
        legend.text.align=0)

legend <- cowplot::get_legend(ggplot(legend_data, aes(x=Superfamily, y=Percentage, fill=Treatment)) + geom_bar(stat="identity") +
                                          scale_fill_manual(breaks=levels(legend_data$Treatment),
                                                            labels=c(expression(paste("C", H[4], ":AMS")),
                                                                     expression(paste("C", H[3], "OH:AMS")),
                                                                     expression(paste("C", H[4], ":NMS")),
                                                                     expression(paste("C", H[3], "OH:NMS"))),
                                                            values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6"),
                                                            guide=guide_legend(ncol=2, title.position="top")) +
                                          theme(legend.position="bottom",
                                                legend.title=element_text(size=9, color="black", hjust=0.5),
                                                legend.text=element_text(size=8, color="black"),
                                                legend.text.align=0,
                                                legend.key.size=unit(1,"line")))

# Assemble ####
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, legend,
             layout_matrix=rbind(c(1,2,3),
                                 c(4,5,6),
                                 c(NA,7,NA)),
             ncol=3, nrow=3, widths=c(1,0.8,1), heights=c(0.4,0.4,0.1))

Fig3 <- arrangeGrob(plot1, plot2, plot3, plot4, plot5, plot6, legend,
                         layout_matrix=rbind(c(1,2,3),
                                             c(4,5,6),
                                             c(NA,7,NA)),
                         ncol=3, nrow=3, widths=c(1,0.8,1), heights=c(0.4,0.4,0.1))

ggsave("Fig3_CH4vCH3OH_2axis.jpg", Fig3, width=6.25, height=5.25, units="in", dpi=300)

rm(plot1, plot2, plot3, plot4, plot5, plot6, Fig3, legend_data)

#### FIGURE 4: CREATE A LEGEND FOR CYTOSCAPE PLOTS ####
fake_data <- data.frame(
  Treatment=c(rep("CH4",50),rep("MeOH",50)),
  Sample=c(1:100),
  log2FC=runif(100, min=-3, max=3)
)
fake_data$Sample <- as.factor(fake_data$Sample)

ggplot(fake_data, aes(x=Sample,y=Treatment)) + geom_tile(aes(fill=Value)) +
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  guides(fill=guide_colorbar(title.position="top",
                             title.hjust=0.5)) +
  theme(legend.position="bottom",
        legend.title=element_text(size=8),
        legend.text=element_text(size=8))

heat_legend <- cowplot::get_legend(ggplot(fake_data, aes(x=Sample,y=Treatment)) + geom_tile(aes(fill=Value)) +
                                     scale_fill_gradient2(low="red", mid="white", high="blue") +
                                     guides(fill=guide_colorbar(title.position="top",
                                                                title.hjust=0.5)) +
                                     theme(legend.position="bottom",
                                           legend.title=element_text(size=8)))

ggsave("~/heat_legend.jpg", heat_legend, dpi=300, height=1, width=1.5, units="in")

rm(fake_data, heat_legend)



#### FIGURE 5: AMMONIUM VS NITRATE ####
# (a, d) STACKED VOLCANO ####
plot7a <- ggplot() +
  geom_point(data=rnaseq_plot_data$AMSvNMS_Methane, aes(x=logFC, y=log10q),
             color=rnaseq_plot_data[["AMSvNMS_Methane"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(y="-log10 q value", tag="a") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,55), expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9, color="white"))

plot7b <- ggplot() +
  geom_point(data=metab_plot_data$AMSvNMS_Methane, aes(x=logFC, y=log10q),
             color=metab_plot_data[["AMSvNMS_Methane"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="a") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", color="white"),
        axis.text=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="white"))

plot7a <- ggplotGrob(plot7a)
plot7b <- ggplotGrob(plot7b)
maxWidth = grid::unit.pmax(plot7a$widths[2:5], plot7b$widths[2:5])
plot7a$widths[2:5] <- as.list(maxWidth)
plot7b$widths[2:5] <- as.list(maxWidth)

plot7 <- arrangeGrob(plot7a, plot7b, nrow=2, heights=c(0.45,0.55))
rm(plot7a, plot7b)

plot10a <- ggplot() +
  geom_point(data=rnaseq_plot_data$AMSvNMS_Methanol, aes(x=logFC, y=log10q),
             color=rnaseq_plot_data[["AMSvNMS_Methanol"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(y="-log10 q value", tag="d") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,55), expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=9, color="white"))

plot10b <- ggplot() +
  geom_point(data=metab_plot_data$AMSvNMS_Methanol, aes(x=logFC, y=log10q),
             color=metab_plot_data[["AMSvNMS_Methanol"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="d") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", color="white"),
        axis.text=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="white"))

plot10a <- ggplotGrob(plot10a)
plot10b <- ggplotGrob(plot10b)
maxWidth = grid::unit.pmax(plot10a$widths[2:5], plot10b$widths[2:5])
plot10a$widths[2:5] <- as.list(maxWidth)
plot10b$widths[2:5] <- as.list(maxWidth)

plot10 <- arrangeGrob(plot10a, plot10b, nrow=2, heights=c(0.45,0.55))
rm(plot10a, plot10b)

# (b) RNA COGS - AMS VS NMS - METHANE ####
plot_data <- rnaseq_numbers$AMSvNMS_Methane
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("NMS", "Fill", "AMS"))

plot8 <- ggplot() +
  geom_bar(data=plot_data, aes(x=COG, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$COG))) +
  scale_y_continuous(limits=c(0,25), breaks=c(0,6.5,12.5,18.5,25), labels=c(25,18.5,12.5,6.5,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,6.5,12.5,18.5,25))) +
  theme_bw() + plot_theme +
  labs(y="% of DEGs", x="COG", tag="b") +
  scale_fill_manual(values=c("#70AD47", "white", "#ED7D31"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkorange3"),
        axis.ticks.x.top=element_line(color="darkorange3"),
        axis.line.x.top=element_line(color="darkorange3"),
        axis.text.x.bottom=element_text(color="darkgreen"),
        axis.ticks.x.bottom=element_line(color="darkgreen"),
        axis.line.x.bottom=element_line(color="darkgreen"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (e) RNA COGS - AMS vs NMS - METHANOL ####
plot_data <- rnaseq_numbers$AMSvNMS_Methanol
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("NMS", "Fill", "AMS"))

plot11 <- ggplot() +
  geom_bar(data=plot_data, aes(x=COG, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$COG))) +
  scale_y_continuous(limits=c(0,25), breaks=c(0,6.5,12.5,18.5,25), labels=c(25,18.5,12.5,6.5,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,6.5,12.5,18.5,25))) +
  theme_bw() + plot_theme +
  labs(y="% of DEGs", x="COG", tag="e") +
  scale_fill_manual(values=c("#2e74b6", "white", "#FFC000"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkgoldenrod3"),
        axis.ticks.x.top=element_line(color="darkgoldenrod3"),
        axis.line.x.top=element_line(color="darkgoldenrod3"),
        axis.text.x.bottom=element_text(color="#2e74b6"),
        axis.ticks.x.bottom=element_line(color="#2e74b6"),
        axis.line.x.bottom=element_line(color="#2e74b6"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

# (c) METAB COGS - AMS VS NMS - METHANE ####
plot_data <- metab_numbers$AMSvNMS_Methane
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("NMS", "Fill", "AMS"))

plot9 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="c") +
  scale_fill_manual(values=c("#70AD47", "white", "#ED7D31"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkorange3"),
        axis.ticks.x.top=element_line(color="darkorange3"),
        axis.line.x.top=element_line(color="darkorange3"),
        axis.text.x.bottom=element_text(color="darkgreen"),
        axis.ticks.x.bottom=element_line(color="darkgreen"),
        axis.line.x.bottom=element_line(color="darkgreen"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (f) METAB COGS - AMS VS NMS - METHANOL ####
plot_data <- metab_numbers$AMSvNMS_Methanol
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("NMS", "Fill", "AMS"))

plot12 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="f") +
  scale_fill_manual(values=c("#2e74b6", "white", "#FFC000"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkgoldenrod3"),
        axis.ticks.x.top=element_line(color="darkgoldenrod3"),
        axis.line.x.top=element_line(color="darkgoldenrod3"),
        axis.text.x.bottom=element_text(color="#2e74b6"),
        axis.ticks.x.bottom=element_line(color="#2e74b6"),
        axis.line.x.bottom=element_line(color="#2e74b6"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)


# Assemble ####
grid.arrange(plot7, plot8, plot9, plot10, plot11, plot12, legend,
             layout_matrix=rbind(c(1,2,3),
                                 c(4,5,6),
                                 c(NA,7,NA)),
             ncol=3, nrow=3, widths=c(1,0.8,1), heights=c(0.4,0.4,0.1))


Fig5 <- arrangeGrob(plot7, plot8, plot9, plot10, plot11, plot12,  legend,
                          layout_matrix=rbind(c(1,2,3),
                                              c(4,5,6),
                                              c(NA,7,NA)),
                          ncol=3, nrow=3, widths=c(1,0.8,1), heights=c(0.4,0.4,0.1))

ggsave("Fig5_2axis.jpg", Fig5, width=6.25, height=5.25, units="in", dpi=300)

rm(Fig5, plot7, plot8, plot9, plot10, plot11, plot12, legend, maxWidth)


#### FIGURE S1: UNDETECTED PLOT ####
temp <- rnaseq_counts$tpm
temp[temp > 0] <- 1

undetect_plot_1 <- colSums(temp)
undetect_plot_1 <- 3794-undetect_plot_1

undetect_plot_1 <- data.frame(
  Sample=rnaseq_sample_data$SampleID,
  Treatment=rnaseq_sample_data$Treatment,
  Counts=undetect_plot_1)

temp <- metab_counts$raw
temp[temp > 0] <- 1

undetect_plot_2 <- colSums(temp, na.rm=TRUE)
undetect_plot_2 <- 341-undetect_plot_2

undetect_plot_2 <- data.frame(
  Sample=metab_sample_data$SampleID,
  Treatment=metab_sample_data$Treatment,
  Counts=undetect_plot_2)

undetect_plot_1 <- ggplot(undetect_plot_1, aes(x=Treatment, y=Counts)) + geom_boxplot() +
  labs(x="Treatment", y="# of Undetected Transcripts", tag="a") +
  theme_bw() + plot_theme +
  scale_x_discrete(labels=c(expression(paste("C", H[4], ":AMS")),
                            expression(paste("C", H[3], "OH:AMS")),
                            expression(paste("C", H[4], ":NMS")),
                            expression(paste("C", H[3], "OH:NMS")))) +
  scale_y_continuous(limits=c(0, 55), expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  
undetect_plot_2 <- ggplot(undetect_plot_2, aes(x=Treatment, y=Counts)) + geom_boxplot() +
  labs(x="Treatment", y="# of Undetected Metabolites", tag="b") +
  theme_bw() + plot_theme +
  scale_x_discrete(labels=c(expression(paste("C", H[4], ":AMS")),
                            expression(paste("C", H[3], "OH:AMS")),
                            expression(paste("C", H[4], ":NMS")),
                            expression(paste("C", H[3], "OH:NMS")))) +
  scale_y_continuous(limits=c(0, 120), expand=c(0,0)) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

grid.arrange(undetect_plot_1, undetect_plot_2, nrow=1)

undetected_plot <- arrangeGrob(undetect_plot_1, undetect_plot_2, nrow=1)

ggsave("FigS1_detection_plot.jpg", undetected_plot, dpi=300, width=5, height=4, units="in")

rm(undetect_plot_1, undetect_plot_2, undetected_plot, temp)

#### FIGURE S2: SPLS-DA ANALYSIS ####
# (a) Transcriptome 
plot_data <- as.data.frame(rnaseq_splsda[["model"]][["variates"]][["X"]])
plot_data$Treatment <- rnaseq_sample_data$Treatment

# from https://stackoverflow.com/questions/24268843/ggplot2-stat-ellipse-draw-ellipses-around-multiple-groups-of-points
library(ellipse)
centroids <- aggregate(cbind(comp1,comp2)~Treatment,plot_data,mean)
conf.rgn  <- do.call(rbind,lapply(unique(plot_data$Treatment),function(t)
  data.frame(Treatment=as.character(t),
             ellipse(cov(plot_data[plot_data$Treatment==t,1:2]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))

rnaseq_splsda[["model"]]$explained_variance

comp12 <- ggplot(plot_data, aes(x=comp1, y=comp2, color=Treatment, group=Treatment)) + geom_point(size=3) +
  geom_path(data=conf.rgn) +
  theme_bw() + plot_theme +
  scale_color_manual(breaks=levels(plot_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme(legend.position="none") +
  labs(x="Component 1 (35.6%)", y="Component 2 (14.0%)", tag="a", title="Transcriptome")

centroids <- aggregate(cbind(comp1,comp3)~Treatment,plot_data,mean)
conf.rgn  <- do.call(rbind,lapply(unique(plot_data$Treatment),function(t)
  data.frame(Treatment=as.character(t),
             ellipse(cov(plot_data[plot_data$Treatment==t,c(1,3)]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))

comp13 <- ggplot(plot_data, aes(x=comp1, y=comp3, color=Treatment, group=Treatment)) + geom_point(size=3) +
  geom_path(data=conf.rgn) +
  theme_bw() + plot_theme +
  scale_color_manual(breaks=levels(plot_data$Treatment),
                   labels=c(expression(paste("C", H[4], ":AMS")),
                            expression(paste("C", H[3], "OH:AMS")),
                            expression(paste("C", H[4], ":NMS")),
                            expression(paste("C", H[3], "OH:NMS"))),
                   values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme(legend.position="none",
        plot.tag=element_text(color="white"),
        plot.title=element_text(color="white")) +
  labs(x="Component 1 (35.6%)", y="Component 3 (12.4%)", tag="x", title="Transcriptome")

rnaseq_splsda[["comp12"]] <- comp12
rnaseq_splsda[["comp13"]] <- comp13

grid.arrange(rnaseq_splsda[["comp12"]], rnaseq_splsda[["comp13"]], nrow=2)

rm(plot_data, comp12, comp13, centroids, conf.rgn)

# (b) Metabolome
plot_data <- as.data.frame(metab_splsda[["model"]][["variates"]][["X"]])
plot_data$Treatment <- metab_sample_data$Treatment

metab_splsda[["model"]]$explained_variance

comp12 <- ggplot(plot_data, aes(x=comp1, y=comp2, color=Treatment, group=Treatment)) + geom_point(size=3) +
  stat_ellipse(aes(color=Treatment), size=1, level = 0.75, alpha=0, type = "norm", geom = "polygon") +
  theme_bw() + plot_theme +
  scale_color_manual(breaks=levels(plot_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x="Component 1 (57.8%)", y="Component 2 (16.9%)", tag="b", title="Metabolome") +
  theme(legend.position="none")

comp13 <- ggplot(plot_data, aes(x=comp1, y=comp3, color=Treatment, group=Treatment)) + geom_point(size=3) +
  stat_ellipse(aes(color=Treatment), size=1, level = 0.75, alpha=0, type = "norm", geom = "polygon") +
  theme_bw() + plot_theme +
  scale_color_manual(breaks=levels(plot_data$Treatment),
                     labels=c(expression(paste("C", H[4], ":AMS")),
                              expression(paste("C", H[3], "OH:AMS")),
                              expression(paste("C", H[4], ":NMS")),
                              expression(paste("C", H[3], "OH:NMS"))),
                     values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  theme(plot.tag=element_text(color="white"),
        legend.position="none",
        plot.title=element_text(color="white")) +
  labs(x="Component 1 (57.8%)", y="Component 3 (5.8%)", tag="x", title="Metabolome")

metab_splsda[["comp12"]] <- comp12
metab_splsda[["comp13"]] <- comp13

grid.arrange(metab_splsda[["comp12"]], metab_splsda[["comp13"]], nrow=2)

legend <- cowplot::get_legend(ggplot(plot_data, aes(x=comp1, y=comp3, color=Treatment, group=Treatment)) +
                                geom_point() +
                                stat_ellipse(aes(color=Treatment),
                                             size=1, level = 0.75, alpha=0, type = "norm", geom = "polygon") +
                                theme_bw() + plot_theme +
                                scale_color_manual(breaks=levels(plot_data$Treatment),
                                                   labels=c(expression(paste("C", H[4], ":AMS")),
                                                            expression(paste("C", H[3], "OH:AMS")),
                                                            expression(paste("C", H[4], ":NMS")),
                                                            expression(paste("C", H[3], "OH:NMS"))),
                                                   values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
                                guides(color=guide_legend(ncol=2, title.position="top")))

grid.arrange(rnaseq_splsda[["comp12"]], metab_splsda[["comp12"]],
                           rnaseq_splsda[["comp13"]], metab_splsda[["comp13"]],
                           legend,
                           layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
                           heights=c(0.4,0.4,0.15))

splsda.plot <- arrangeGrob(rnaseq_splsda[["comp12"]],
                           metab_splsda[["comp12"]],
                           rnaseq_splsda[["comp13"]],
                           metab_splsda[["comp13"]],
             legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
             heights=c(0.4,0.4,0.15))

ggsave("BG8_splsda_plot.jpg",
       splsda.plot, width=6.5, height=6, units="in", dpi=300) 


splsda1 <- ggplotGrob(rnaseq_splsda[["comp12"]])
splsda2 <- ggplotGrob(rnaseq_splsda[["comp13"]])
maxWidth = grid::unit.pmax(splsda1$widths[2:5], splsda2$widths[2:5])
splsda1$widths[2:5] <- as.list(maxWidth)
splsda2$widths[2:5] <- as.list(maxWidth)

splsda3 <- ggplotGrob(metab_splsda[["comp12"]])
splsda4 <- ggplotGrob(metab_splsda[["comp13"]])
maxWidth = grid::unit.pmax(splsda3$widths[2:5], splsda4$widths[2:5])
splsda3$widths[2:5] <- as.list(maxWidth)
splsda4$widths[2:5] <- as.list(maxWidth)

grid.arrange(splsda1, splsda3,
             splsda2, splsda4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

splsda.plot <- arrangeGrob(splsda1, splsda3,
                           splsda2, splsda4, legend,
                           layout_matrix=rbind(c(1,2),c(3,4),c(5,5)), heights=c(0.4,0.4,0.15))

ggsave("~/BG8_splsda_plot.jpg",
       splsda.plot, width=5, height=5.5, units="in", dpi=300)

rm(class, data, list.keepX, select.keepX, splsda.srbct, plot_data, error, ncomp, legend
   tune.splsda.srbct, comp12, comp13, splsda1, splsda2, splsda3, splsda4, maxWidth, splsda.plot)



########### SUPPLEMENTAL BAR CHARTS: PREPARE DATA FOR PLOTTING ####
plot_key <- read.csv("raw_data/supplemental_plots_info.csv")

# Add mean and sd TPKM values to the rnaseq data.
plot_key_rna <- subset(plot_key, is.na(MetaboliteID)=="TRUE")
plot_key_rna_avg <- merge(plot_key_rna, rnaseq_counts$avg, by="Accession", all=FALSE)
plot_key_rna_sd <- merge(plot_key_rna, rnaseq_counts$sds, by="Accession", all=FALSE)

plot_key_rna_avg <- melt(plot_key_rna_avg)
plot_key_rna_sd <- melt(plot_key_rna_sd)

plot_key_rna_avg$sd <- plot_key_rna_sd$value 
plot_key_rna <- plot_key_rna_avg
rm(plot_key_rna_sd, plot_key_rna_avg)

colnames(plot_key_rna) <- c("Accession", "UniProtKB", "MetaboliteID", "Group", "Name", "Treatment", "Mean", "SD")

# Add mean and sd metabolite abundances to the metabolite data.
plot_key_met <- subset(plot_key, is.na(MetaboliteID)=="FALSE")
plot_key_met_avg <- merge(plot_key_met, metab_counts$avg, by="MetaboliteID", all=FALSE)
plot_key_met_sd <- merge(plot_key_met, metab_counts$sds, by="MetaboliteID", all=FALSE)

plot_key_met_avg <- melt(plot_key_met_avg)
plot_key_met_sd <- melt(plot_key_met_sd)

plot_key_met_avg$sd <- plot_key_met_sd$value 
plot_key_met <- plot_key_met_avg
rm(plot_key_met_sd, plot_key_met_avg)

colnames(plot_key_met) <- c("MetaboliteID", "Accession", "UniProtKB", "Group", "Name", "Treatment", "Mean", "SD")
plot_key_met <- plot_key_met[,c("Accession", "UniProtKB", "MetaboliteID", "Group", "Name", "Treatment", "Mean", "SD")]

# Put them back together into a single file for plotting.
plot_key <- rbind(plot_key_rna, plot_key_met)
rm(plot_key_rna, plot_key_met)
levels(plot_key$Treatment)

# Remove error bars that will go below zero.
plot_key$SD[plot_key$SD > plot_key$Mean] <- plot_key$Mean[plot_key$SD > plot_key$Mean]

#### FIGURE S3: STRESS RESPONSE ####
# General stress
plot1 <- ggplot(subset(plot_key, Group=="General"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance (TPM)", title="csrA and sigma factors") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot1 <- arrangeGrob(plot1, blank, blank, ncol=3, widths=c(4,1,1))

# Hopanoid biosynthesis  
plot2 <- ggplot(subset(plot_key, Group=="Hopanoid biosynthesis"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance (TPM)", title="Hopanoid biosynthesis") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

# Oxidative stress
plot3 <- ggplot(subset(plot_key, Group=="Oxidative stress"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance (TPM)", title="Oxidative stress") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

grid.arrange(plot1, plot2, plot3, nrow=3, heights=c(1,1,2))

Fig_S3 <- arrangeGrob(plot1, plot2, plot3, nrow=3, heights=c(1,1,2))
ggsave("~/Figure_S3_stress.jpg", Fig_S3, width=8, height=6, units="in", dpi=300)

#### FIGURE S4: METAL TRANSPORT, NON-RIBOSOMAL PEPETIDE SYNTHASES ####
plot4 <- ggplot(subset(plot_key, Group=="Metal transport"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Metal transporters") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot5 <- ggplot(subset(plot_key, Group=="Peptide synthase"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Non-ribosomal peptide synthases") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot5 <- arrangeGrob(plot5, blank, blank, blank, ncol=4, widths=c(3.2,1,1,1))

grid.arrange(plot4, plot5, nrow=2, heights=c(1,0.5))

Fig_S4 <- arrangeGrob(plot4, plot5, nrow=2, heights=c(1,0.5))
ggsave("~/Figure_S4_metal_NRPS.jpg", Fig_S4, width=8, height=5, units="in", dpi=300)

#### FIGURE S5: OXIDATIVE PHOSPHORYLATION ####
plot6 <- ggplot(subset(plot_key, Group=="Cytochromes"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Cytochromes") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot7 <- ggplot(subset(plot_key, Group=="NADH oxidoreductase"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="NADH oxidoreductase") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

grid.arrange(plot6, plot7, nrow=2, heights=c(1,3))

Fig_S5 <- arrangeGrob(plot6, plot7, nrow=2, heights=c(1,3))
ggsave("~/Figure_S5_oxido_phosphoro.jpg", Fig_S5, width=8.5, height=6.5, units="in", dpi=300)

#### FIGURE S6: GAMMA-GLUTMAYL AMINO ACIDS ####
plot8 <- ggplot(subset(plot_key, Group=="Gamma-glutamyl Aas"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Gamma-glutamyl amino acids") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

ggsave("~/Figure_S6_gammaglutamyl.jpg", plot8, width=8.5, height=4.5, units="in", dpi=300)

#### FIGURE S7: SELECT ITEMS UPREGULATED IN METHANOL ####
# Phospholipids and fatty acids - 7 items
plot14 <- ggplot(subset(plot_key, Group=="Phospholipids" | Group=="Fatty acids"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=7) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Phospholipids and fatty acids") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

# Aromatic amino acids - 10 items
plot15 <- ggplot(subset(plot_key, Group=="Aromatic AA"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Aromatic amino acids") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

# Branched chain amino acids - 11 items
plot16 <- ggplot(subset(plot_key, Group=="BCAA"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Branched chain amino acids") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot17 <- ggplot(subset(plot_key, Group=="Other methanol"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Histidine betaine and histidinol") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot17 <- arrangeGrob(plot17, blank, ncol=2, widths=c(2.1,4))

grid.arrange(plot14, plot15, plot16, plot17, nrow=4, heights=c(1,1,1,1))
Fig_S <- arrangeGrob(plot14, plot15, plot16, plot17, nrow=4, heights=c(1,1,1,1))

ggsave("~/Figure_S_methanol_metabs.jpg", Fig_S, width=8.5, height=6, units="in", dpi=300)

#### FIGURE S8, S9: GLYCOLYTIC PATHWAYS ####
# EMP pathway - 14 boxes
plot9 <- ggplot(subset(plot_key, Group=="EMP"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=5) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Embden-Meyerhof-Parnas pathway") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

# Shared EDD/PPP pathway - 4 boxes
plot10 <- ggplot(subset(plot_key, Group=="PPP/EDD"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=5) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Shared enzymes (PPP and EDD)") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot10 <- arrangeGrob(plot10, blank, ncol=2, widths=c(4,1))

# PPP pathway - 10 boxes
plot11 <- ggplot(subset(plot_key, Group=="PPP"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=5) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Pentose phosphate pathway") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

# EDD pathway - 4 boxes
plot12 <- ggplot(subset(plot_key, Group=="EDD"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=5) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Enter-Doudoroff pathway") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

plot12 <- arrangeGrob(plot12, blank, ncol=2, widths=c(4,1))

plot9
grid.arrange(plot10, plot11, plot12, nrow=3, heights=c(1,2,1))

ggsave("~/Figure_S8_EMP_pathway.jpg", plot9, width=8, height=5, units="in", dpi=300)

Fig_S7 <- arrangeGrob(plot10, plot11, plot12, nrow=3, heights=c(1,2,1))
ggsave("~/Figure_S7_PPP.EDD_pathways.jpg", Fig_S7, width=8, height=6, units="in", dpi=300)


#### FIGURE S10: CARBON ENZYMES ####
plot18 <- ggplot(subset(plot_key, Group=="Carbon"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Methane cycle enzymes") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

ggsave("~/Fig_S_methane_enzymes.jpg", plot18, width=8.5, height=4.5, units="in", dpi=300)

#### FIGURE S12: NITROGEN METABOLISM ####
plot13 <- ggplot(subset(plot_key, Group=="Nitrogen metabolism"), aes(x=Treatment, y=Mean, fill=Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
  theme_bw() + 
  scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
  labs(x=NULL, y="Abundance", title="Nitrogen metabolism") +
  theme(axis.text.x=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        plot.title=element_text(size=11, face="bold"),
        strip.text=element_text(size=8, color="black"),
        legend.text=element_text(hjust=0),
        legend.position="none")

ggsave("~/Figure_S_Nitrogen_metabolism.jpg", plot13, width=8.5, height=4, units="in", dpi=300)

#### CREATE LEGENDS FOR SUPPLEMENTAL BAR CHARTS ####
legend2col <- cowplot::get_legend(ggplot(subset(plot_key, Group=="Carbon"), aes(x=Treatment, y=Mean, fill=Treatment)) +
                                    geom_bar(stat = "identity", position = "dodge") +
                                    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
                                    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                                    facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
                                    theme_bw() + 
                                    scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                                                      labels=c(expression(paste("C", H[4], ":AMS")),
                                                               expression(paste("C", H[3], "OH:AMS")),
                                                               expression(paste("C", H[4], ":NMS")),
                                                               expression(paste("C", H[3], "OH:NMS"))),
                                                      values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
                                    guides(fill=guide_legend(ncol=2)) +
                                    labs(x=NULL, y="Abundance", title="Methane cycle enzymes") +
                                    theme(axis.text.x=element_blank(),
                                          panel.grid=element_blank(),
                                          axis.text.y=element_text(size=8, color="black"),
                                          plot.title=element_text(size=11, face="bold"),
                                          strip.text=element_text(size=8, color="black"),
                                          legend.text=element_text(hjust=0),
                                          legend.title=element_text(hjust=0.5, size=10, color="black"),
                                          legend.position="right"))

legend1col <- cowplot::get_legend(ggplot(subset(plot_key, Group=="Carbon"), aes(x=Treatment, y=Mean, fill=Treatment)) +
                                    geom_bar(stat = "identity", position = "dodge") +
                                    geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=0.3) +
                                    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
                                    facet_wrap(~Name, scales="free", labeller = label_wrap_gen(15), ncol=6) +
                                    theme_bw() + 
                                    scale_fill_manual(breaks=levels(rnaseq_sample_data$Treatment),
                                                      labels=c(expression(paste("C", H[4], ":AMS")),
                                                               expression(paste("C", H[3], "OH:AMS")),
                                                               expression(paste("C", H[4], ":NMS")),
                                                               expression(paste("C", H[3], "OH:NMS"))),
                                                      values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6")) +
                                    guides(fill=guide_legend(ncol=1)) +
                                    labs(x=NULL, y="Abundance", title="Methane cycle enzymes") +
                                    theme(axis.text.x=element_blank(),
                                          panel.grid=element_blank(),
                                          axis.text.y=element_text(size=8, color="black"),
                                          plot.title=element_text(size=11, face="bold"),
                                          strip.text=element_text(size=8, color="black"),
                                          legend.text=element_text(hjust=0),
                                          legend.title=element_text(hjust=0.5, size=10, color="black"),
                                          legend.position="right"))

ggsave("~/legend2col.jpg", legend2col, width=2.5, height=1.2, units="in", dpi=300)
ggsave("~/legend1col.jpg", legend1col, width=1.25, height=2, units="in", dpi=300)

########### NORMALIZED DATA PLOT #########
#### FIGURE 3: METHANE VS METHANOL ####
# (a, c) STACKED VOLCANO ####
# Create column identifying significantly differentially abundant items in red.
norm_plot_data <- norm_diff_abund
for(i in c(1:4)){
  norm_plot_data[[i]]$color <- rep("grey")
  norm_plot_data[[i]] <- subset(norm_plot_data[[i]], is.na(logFC)=="FALSE")
  for(p in 1:nrow(norm_plot_data[[i]])){
    if(norm_plot_data[[i]][p,"p.adj"] < 0.01 & abs(norm_plot_data[[i]][p,"logFC"]) > 1){
      norm_plot_data[[i]][p,"color"] <- "red"
    }
  }
}

# Metabolome: Methane vs. methanol, AMS
plot1 <- ggplot() +
  geom_point(data=norm_plot_data$CH4vCH3OH_AMS, aes(x=logFC, y=log10q),
             color=norm_plot_data[["CH4vCH3OH_AMS"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="a") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,8), expand=c(0,0)) +
  theme(legend.position="none",
        plot.tag=element_text(face="bold", color="black"),
        axis.text=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))

# Metabolome: Ammonium vs nitrate, in methanol
plot3 <- ggplot() +
  geom_point(data=norm_plot_data$AMSvNMS_Methanol, aes(x=logFC, y=log10q),
             color=norm_plot_data[["AMSvNMS_Methanol"]]$color, size=1) +
  theme_bw() + plot_theme +
  labs(x="log2 Fold Change", y="-log10 q value", tag="c") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0, color="grey") +
  scale_x_continuous(limits=c(-10,10)) +
  scale_y_continuous(limits=c(0,10), expand=c(0,0)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=8, color="black"),
        axis.title.x=element_text(size=9, color="black"),
        axis.title.y=element_text(size=9, color="black"))


# (c) NORM COGS - CH4 v CH3OH - AMS ####
plot_data <- norm_numbers$CH4vCH3OH_AMS
plot_data$Treatment <- factor(plot_data$Treatment,
                              levels=c("CH3OH", "Fill", "CH4"))

plot2 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="b") +
  scale_fill_manual(values=c("#FFC000", "white", "#ED7D31"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="orangered3"),
        axis.ticks.x.top=element_line(color="orangered3"),
        axis.line.x.top=element_line(color="orangered3"),
        axis.text.x.bottom=element_text(color="darkgoldenrod3"),
        axis.ticks.x.bottom=element_line(color="darkgoldenrod3"),
        axis.line.x.bottom=element_line(color="darkgoldenrod3"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())

rm(plot_data)

# (f) NORM COGS - AMS v NMS - METHANOL ####
plot_data <- norm_numbers$AMSvNMS_Methanol
plot_data$Treatment <- factor(plot_data$Treatment, levels=c("NMS", "Fill", "AMS"))

plot4 <- ggplot() +
  geom_bar(data=plot_data, aes(x=Superfamily, y=Percentage, fill=Treatment), stat="identity") +
  coord_flip() +
  scale_x_discrete(limits=rev(levels(plot_data$Superfamily))) +
  scale_y_continuous(limits=c(0,100), breaks=c(0,25,50,75,100), labels=c(100,75,50,25,0), expand=c(0,0),
                     sec.axis = dup_axis(labels=c(0,25,50,75,100))) +
  theme_bw() + plot_theme +
  labs(y="% of DAMs", x="Superfamily", tag="d") +
  scale_fill_manual(values=c("#2e74b6", "white", "#FFC000"), guide=guide_legend(ncol=1, title.position="top")) +
  theme(legend.position="none",
        panel.background = element_rect(colour = "black", size=0.5, fill=NA),
        axis.text.x.top=element_text(color="darkgoldenrod3"),
        axis.ticks.x.top=element_line(color="darkgoldenrod3"),
        axis.line.x.top=element_line(color="darkgoldenrod3"),
        axis.text.x.bottom=element_text(color="#2e74b6"),
        axis.ticks.x.bottom=element_line(color="#2e74b6"),
        axis.line.x.bottom=element_line(color="#2e74b6"),
        axis.text.y=element_text(size=7),
        axis.title.x.top=element_blank())


rm(plot_data)

# Create a legend ####
legend_data <- rbind(norm_numbers$AMSvNMS_Methane, norm_numbers$CH4vCH3OH_AMS)
legend_data <- subset(legend_data, Treatment !="Fill")
legend_data$Treatment <- factor(legend_data$Treatment)

ggplot(legend_data, aes(x=Superfamily, y=Percentage, fill=Treatment)) + geom_bar(stat="identity") +
  scale_fill_manual(breaks=levels(legend_data$Treatment),
                    labels=c(expression(paste("C", H[4], ":AMS")),
                             expression(paste("C", H[3], "OH:AMS")),
                             expression(paste("C", H[4], ":NMS")),
                             expression(paste("C", H[3], "OH:NMS"))),
                    values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6"),
                    guide=guide_legend(ncol=2, title.position="top")) +
  theme(legend.position="bottom",
        legend.title=element_text(size=10, color="black", hjust=0.5),
        legend.text=element_text(size=9, color="black"),
        legend.text.align=0)

legend <- cowplot::get_legend(ggplot(legend_data, aes(x=Superfamily, y=Percentage, fill=Treatment)) + geom_bar(stat="identity") +
                                scale_fill_manual(breaks=levels(legend_data$Treatment),
                                                  labels=c(expression(paste("C", H[4], ":AMS")),
                                                           expression(paste("C", H[3], "OH:AMS")),
                                                           expression(paste("C", H[4], ":NMS")),
                                                           expression(paste("C", H[3], "OH:NMS"))),
                                                  values=c("#ED7D31", "#FFC000", "#70AD47", "#2e74b6"),
                                                  guide=guide_legend(ncol=2, title.position="top")) +
                                theme(legend.position="bottom",
                                      legend.title=element_text(size=9, color="black", hjust=0.5),
                                      legend.text=element_text(size=8, color="black"),
                                      legend.text.align=0,
                                      legend.key.size=unit(1,"line")))

plot1 <- ggplotGrob(plot1)
plot3 <- ggplotGrob(plot3)
plot1$widths <- plot3$widths

# Assemble ####
grid.arrange(plot1, plot2, plot3, plot4, legend,
             layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
             widths=c(1, 1), heights=c(0.4,0.4,0.1))

supp_panels <- arrangeGrob(plot1, plot2, plot3, plot4, legend,
                           layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),
                           widths=c(1, 1), heights=c(0.4,0.4,0.1))

ggsave("FigS14_relativized_metabs.jpg", supp_panels, width=5, height=5.25, units="in", dpi=300)

rm(plot1, plot2, plot3, plot4, supp_panels, legend, legend_data)
