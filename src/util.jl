minimum_nonmissing(x) = all(ismissing.(x)) ? missing : minimum(skipmissing(x))
minimum_nonmissing(x...) = minimum_nonmissing(x)