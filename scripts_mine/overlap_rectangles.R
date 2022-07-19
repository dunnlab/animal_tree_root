
library(ggplot2)
require(gridExtra)


df = overlap_rects

#xmin_l <- df$xmin
#xmax_l <- df$xmax
#ymin_l <- df$ymin
#ymax_l <- df$ymax


#plot_all <- ggplot(df, aes_string(xmin = df$xmin, xmax =df$xmax, ymin = df$ymin, ymax =df$ymax )) +
#	geom_rect(color = NA, fill = "pink", alpha = 0.4) #+ theme_void() + theme(plot.background = element_rect(fill = "black"))

#df_int <- df[c(4,5),]

#plot_intersecting <- ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
#	geom_rect(fill = "pink", alpha = 0.4) #+ theme_void() + theme(plot.background = element_rect(fill = "black"))



p_list <- list()
i=0

for (val in seq(1,nrow(df),2))
{
	df_int <- df[c(val,val+1),]
	title=paste(df_int$matrix_1[1],df_int$matrix_2[1])
	i = i+1
	#nam <- paste("plot", val, sep = "")
	
	#assign(nam, ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
				 	#geom_rect(fill = "pink", alpha = 0.4) )#+ theme_void() + theme(plot.background = element_rect(fill = "black"))
	p_list[[i]] <- ggplot(df_int, aes_string(xmin = df_int$xmin, xmax =df_int$xmax, ymin = df_int$ymin, ymax =df_int$ymax )) +
		geom_rect(fill = c(1,2), alpha = 0.4) +
		  ggtitle(title) + theme(plot.title = element_text(size = 4))
}

do.call(grid.arrange,p_list)




#plot_em = ggplot(df)+geom_blank()

#grid.arrange(plot2, plot_em, plot4, plot7, ncol=floor(sqrt(nrow(df))), nrow=ceiling(sqrt(nrow(df))))

#lapply(dflist, function(df) {
	# Do some complex operations on each data frame, df
	# More steps
	
	# Make sure the last thing is NULL. The last statement within the function will be
	# returned to lapply, which will try to combine these as a list across all data frames.
	# You don't actually care about this, you just want to run the function.
	
	
#	NULL
#})

#grid.arrange(plot2, plot_em, plot4, plot7, ncol=2, nrow=2)

#plot + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))