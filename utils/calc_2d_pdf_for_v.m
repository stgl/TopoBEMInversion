function pdf = calc_2d_pdf_for_v(X, vs, vn, sig_vs, sig_vn, cov_v)

pdf = zeros(length(X(:,1)),length(X(1,:)));

n = length(vs);
for(i=1:n)
  pdf = pdf + mvnpdf(X,[vs(i) vn(i)],[sig_vs cov_v;cov_v sig_vn]);
end

pdf = pdf ./ n;
