% Define the filename
filename = 'ch1_MVC_envelope.csv';

% Combine the processed data into a single matrix
% Each column represents a variable (adjust this order as needed)
data_to_save = [sEMG_envelope1];

% Save the matrix to a CSV file with headers
header = {'ch1_MVC_envelope'};  % Define headers

% Write headers and data to the CSV file
writecell(header, filename);              % Write the header
writematrix(data_to_save, filename, 'WriteMode', 'append');  % Append data below the header
