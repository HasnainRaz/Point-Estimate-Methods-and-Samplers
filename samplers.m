function [histobj, lower, higher] = samplers(name, x, fx, iterations, parameters)
    switch name
        case 'mh'
            if nargin < 5
                parameters = [0, 1];
            end
            init_x = 0.0001;
            accepted = zeros(1, iterations);
            counter = 0;
            alpha = -1;
            beta = 1;
            prob_proposal = 0;
            prob_current = 0;
            for i = 1:1:iterations
                proposal = init_x + alpha + (beta-alpha)*rand(1);
                if isa(fx, 'char')
                prob_proposal = distri(proposal, fx, parameters, 'pdf');
                prob_current = distri(init_x, fx, parameters, 'pdf');
                else
                prob_proposal = interp1(x, fx, proposal);
                prob_current = interp1(x, fx, init_x);
                end
                acceptance_prob = min([1, prob_proposal/prob_current]);
                chance = rand(1);
                if chance < acceptance_prob
                    init_x = proposal;
                    counter = counter + 1;
                    accepted(counter) = proposal;
                end
            end
            
            accepted = accepted(1:counter);
            lower = min(accepted);
            higher = max(accepted);
            
            histobj = accepted;
            
        case 'lhs'
            segment_size = 1/iterations;
            segment_start = 0;
            sample_probability = zeros(1, iterations);
            cdf = cumtrapz(x, fx);
            non_mono_indices = find(cdf == 0);
            non_mono_indices = non_mono_indices(1:end - 1);
            non_mono_indices = [non_mono_indices, find(diff(cdf) == 0)];
            for i=1:1:iterations
                segment_end = segment_start + segment_size;
                sample_probability(i) = segment_start + (segment_end - segment_start)*rand();
                segment_start = segment_end;
            end
            cdf(non_mono_indices) = [];
            x(non_mono_indices) = [];
            samples = interp1(cdf, x, sample_probability);
            lower = min(samples);
            higher = max(samples);
            histobj = samples;
    end
    
end