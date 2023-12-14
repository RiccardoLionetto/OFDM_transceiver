function output = lfsr_training(output_length)
% A linear feedback shift register (LFSR) which outputs a PN sequence of length output_length.
% The current LFSR has a period length of 255, but the polynomial can easily be changed for a longer one.

polynomial = [1 0 1 1 1 0 0 0]';

% All memories are initialized with ones
state = ones(size(polynomial));

output = zeros(output_length, 1);

for i = 1:output_length
    output(i) = state(1);
    feedback = mod(sum(state .* polynomial), 2);
    state = circshift(state, -1);
    state(end) = feedback;
end
end