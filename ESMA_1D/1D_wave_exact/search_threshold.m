function search_threshold(nu, z)
vareps = 10^-7;

y1 = calculate_Mainardi_function(z, nu, 1000, vareps)
y2 = calculate_Mainardi_asymtotic(z, nu)

end              