library(RWsearch)
crandb_down()
s_crandb("extreme", "value", select="D", mode="and")
       
p_text(s_crandb("extreme", "value", select="D", mode="and"))

s_crandb("peak", "over", "threshold", select="D", mode="and")


p_text(s_crandb("peaks", "over", "threshold", select="D", mode="and"))


pkgtolook <- c(s_crandb("peaks", "over", "threshold", select="D", mode="and"),
               s_crandb("extreme", "value", select="D", mode="and"))

table(pkgtolook)[table(pkgtolook) > 1]
pkgtolook <- sort(unique(pkgtolook))

p_down(pkgtolook, dir="pkg-to-look-at")

p_down(c("QRM", "ReIns"), dir=".")
p_down("ExtremeRisks", dir=".")
p_down("RTDE", dir="pkg-to-look-at")


x <- s_crandb("extreme", "Stan", select="D", mode="and", sensitive=TRUE)


x[!x %in% pkgtolook]

p_text(x[!x %in% pkgtolook])


x <- s_crandb("extreme", "INLA", select="D", mode="and")
