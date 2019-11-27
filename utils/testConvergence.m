function test_results = testConvergence(samples, start_sample_number, number_of_samples)

[n, ~] = size(samples);
[m, ~] = size(samples{1});

mat_samples = reshape(cell2mat(samples),m,n)';

test_results = zeros(m,1);

for(i=1:m)
    sample1 = mat_samples(start_sample_number:start_sample_number+number_of_samples,i);
    sample2 = mat_samples(start_sample_number+number_of_samples+1:start_sample_number+number_of_samples*2, i);
    test_results(i) = kstest2(sample1, sample2);
end
