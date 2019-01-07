
AvgBedgraph = function (df) {

avg = vector(length=dim(df)[1])

# get the average of col_4 ... col_n

cov_values = df[,4:dim(df)[2]]
cov_values = t(cov_values)

for (i in 1:dim(cov_values)[2]) {

        # only retrieve the bed record columns
        avg[i] = mean(cov_values[,i])
}

df$avg = avg

return(df)
}

AvgSalmonQuant = function (df) {

avg = vector(length=dim(df)[1])

# get the average of col_4 ... col_n

cov_values = df[,2:dim(df)[2]]
cov_values = t(cov_values)

for (i in 1:dim(cov_values)[2]) {

        # only retrieve the bed record columns
        avg[i] = mean(cov_values[,i])
}

df$avg = avg

return(df)
}
