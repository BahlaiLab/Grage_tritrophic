#run Konza_analysis.R, Hubbard_analysis.R and SBC_analysis.R first

#three timeseries figure

timeseriesplot.hubbarddata1<-timeseriesplot.hubbarddata+theme(legend.position = "inside",
  legend.position.inside = c(0.1, 0.75))
timeseriesplot.konzadata1<-timeseriesplot.konzadata+theme(legend.position="none")
timeseriesplot.sbcdata1<-timeseriesplot.sbcdata+theme(legend.position="none")

all.timeseries.site<-plot_grid(timeseriesplot.hubbarddata1, timeseriesplot.konzadata1, timeseriesplot.sbcdata1,
                                     ncol=1,  rel_heights=c(1,1,1), labels="AUTO", axis="b", align="v")

all.timeseries.site

pdf("figures/timeseries.site.pdf", width = 8, height = 8)
all.timeseries.site
dev.off()

#three heatmap figure

all_heat.hubbarddata1<-all_heat.hubbarddata+ theme(legend.text = element_blank(),
                                           legend.title = element_blank(),
                                           legend.key = element_blank(),
                                           legend.background = element_rect(fill = "transparent", colour = NA),
                                           axis.title.x  = element_text(color = "transparent"), 
                                           axis.title.x.top  = element_text(color = "black"))+
  guides(color = guide_none(), fill = guide_none(), shape = guide_none(), linetype = guide_none())

all_heat.konzadata1<-all_heat.konzadata+ theme(axis.title.x  = element_text(color = "transparent"),
                                               axis.title.x.top  = element_text(color = "transparent"))

all_heat.sbcdata1<-all_heat.sbcdata+ theme(legend.text = element_blank(),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA), 
      axis.title.x.top  = element_text(color = "transparent"))+
  guides(color = guide_none(), fill = guide_none(), shape = guide_none(), linetype = guide_none())

all_heat.site<-plot_grid(all_heat.hubbarddata1, all_heat.konzadata1, all_heat.sbcdata1,
                               ncol=1,  rel_heights=c(1,1,1), labels="AUTO", axis="b", align="v")

all_heat.site

pdf("figures/all_heat.site.pdf", width = 8, height = 8)
all_heat.site
dev.off()



allccf.sbcdata
allderived.sbcdata

complete_analysis.sbcdata<-plot_grid(timeseriesplot.sbcdata, all_heat.sbcdata, allccf.sbcdata, allderived.sbcdata, 
                                     ncol=2,  rel_heights=c(1,1,1,1), labels="AUTO", axis="b", align="v")

complete_analysis.sbcdata

pdf("figures/sbcfig1.pdf", width = 14.5, height = 9)
complete_analysis.sbcdata
dev.off()