using TableReader, DataFrames

grat = readcsv("data/grat.csv")
rename!(grat, :UNNAMED_1 => :ID)
# convert
y = convert(Matrix{Int64}, grat[!, Not(:ID)])
rs = sum(y; dims = 2)
# Get unique index
unique_ind = findall(x -> !x, nonunique(grat[!, Not(:ID)]))
# unique_ind = findall(x -> !x, nonunique(y, 1))
