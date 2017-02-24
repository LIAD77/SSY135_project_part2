clear; close;

% Set parameters
modulations = [4, 16, 64];
reps = [1, 2, 3];
markers = {'o', 'd', '+'};
colors = {'r', 'b', 'k'};

% Link budget for the various constellations
EbN0_4 = 19.65;
EbN0_16 = 16.64;
EbN0_64 = 14.88;
link_budget = [EbN0_4, EbN0_16, EbN0_64];

for modulation_index = 1:length(modulations)
    modulation_order = modulations(modulation_index);

    figure; hold on;
    for rep_index = 1:length(reps)
        repetitions = reps(rep_index);
        color = colors{rep_index};
        marker = markers{rep_index};

        % Parse the parameters into a filename
        filename = strcat(num2str(modulation_order), '_', num2str(repetitions), '.mat');
        BER_file = strcat('BER', filename);
        EbN0_file = strcat('EbN0', filename);

        % Run the simulations if the results were not found on disk
        if exist(BER_file, 'file') == 0 || exist(EbN0_file, 'file') == 0
            disp(['Simulating ', filename])
            eval('project')
        end

        % Load data from disk
        BER = load(BER_file); BER = BER.BER;
        EbN0 = load(EbN0_file); EbN0 = EbN0.EbN0_sequence;

        % Plot the empirical BER
        plot(EbN0, BER, strcat([color, marker]))

        % Plot the theoretical BER
        BER_theory = zeros(size(EbN0));
        for index = 1:length(EbN0)

            % Convert SNR to watts
            EbN0_W = 10^(EbN0(index) / 10);

            % Gamma distribution otherwise
            b = EbN0_W / repetitions;
            a = EbN0_W / b;
            fun = @(gamma) qamerr(gamma, modulation_order) .* gampdf(gamma, a, b);
            BER_theory(index) = integral(fun, 0, inf);
        end
        plot(EbN0, BER_theory, strcat([color, '-']))
    end

    % Plot the link budget
    budget = link_budget(modulation_index);
    plot([budget, budget], [1, 2e-4], 'LineWidth', 2)

    % Plot the required EbN0 and max EbN0
    BER_req = 2e-4;
    plot([0, 20], [BER_req, BER_req], 'LineWidth', 2)

    legend({'$r=1$ Empirical', '$r=1$ Theory', '$r=2$ Empirical', '$r=2$ Theory', ...
            '$r=3$ Empirical', '$r=3$ Theory', 'Link Budget', 'BER Requirement'}, ...
           'Interpreter', 'LaTex', 'Location', 'SouthWest')

    xlabel('$\frac{E_b}{N_0}$', 'Interpreter', 'LaTex')
    ylabel('BER', 'Interpreter', 'LaTex')
    set(gca, 'FontSize', 15)

    grid on; hold off;
    set(gca,'yscale','log');
end

return
