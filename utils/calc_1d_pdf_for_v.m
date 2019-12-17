function pdf = calc_1d_pdf_for_v(x, mus, sigmas)

n = length(mus);
pdf = zeros(1,length(x));
for(i=1:n)
  pdf = pdf + normpdf(x,mus(i),sigmas(i));
end

pdf = pdf ./ n;
