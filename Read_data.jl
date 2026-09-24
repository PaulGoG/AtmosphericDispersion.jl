# Citirea CSV-urilor cu date tabelate cu care o sa efectuam calculele

cd(@__DIR__)

T_1 = CSV.File("Tabulated_data\\Table_1.csv"; normalizenames=true) |> DataFrame
T_2 = CSV.File("Tabulated_data\\Table_2.csv"; normalizenames=true) |> DataFrame
T_3 = CSV.File("Tabulated_data\\Table_3.csv"; normalizenames=true) |> DataFrame
T_4 = CSV.File("Tabulated_data\\Table_4.csv"; normalizenames=true) |> DataFrame
T_7 = CSV.File("Tabulated_data\\Table_7.csv"; normalizenames=true) |> DataFrame
Cladiri = CSV.File("Tabulated_data\\Buildings.csv"; normalizenames=true) |> DataFrame
Freq = CSV.File("Tabulated_data\\Frequencies.csv"; normalizenames=true) |> DataFrame