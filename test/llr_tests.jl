# raw log file for first replica of the first repeat
fileSp4 = "./test_data/4x20_0_Rep0_out_0.txt"
# David's comparison files
david_fa_Sp4 = "./test_data/4x20_0_Rep0_fa.csv"
david_RM_Sp4 = "./test_data/4x20_0_Rep0_RM.csv"

fa = readdlm(david_fa_Sp4, ',', skipstart = 1)
RM = readdlm(david_RM_Sp4, ',', skipstart = 1)
# only use rows where the replica (the 5th or 7th entry respectively) number is zero
fa = fa[fa[:, 5] .== 0, :]
RM = RM[RM[:, 7] .== 0, :]
RM = sortslices(RM, dims = 1, by = x -> (x[1]))
fa = sortslices(fa, dims = 1, by = x -> (x[8]))

dS0, S0, plaq, a, is_rm, S0_fxa, a_fxa, poly, E_therm, E_meas, E_fxa = parse_llr(fileSp4)

a_fxa_david = Float64[]
S0_fxa_david = Float64[]
for i in Iterators.partition(eachindex(fa[:, 1]), 100)
    append!(a_fxa_david, only(unique(fa[i, 1])))
    append!(S0_fxa_david, only(unique(fa[i, 2])))
end

@testset "Compare LLR parsing with David's code" begin
    @test dS0 == 2only(unique(RM[:, 4]))
    @test S0 == RM[:, 3]
    @test isapprox(plaq, RM[:, 6] ./ RM[:, 5] ./ 6)
    @test a == -RM[:, 2]
    @test abs.(poly) == fa[:, 6]
    @test a_fxa == -a_fxa_david
    @test S0_fxa == S0_fxa_david
end
