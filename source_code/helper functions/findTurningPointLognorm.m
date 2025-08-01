function userValue = findTurningPoint(app, data, ptitle, q)
% inputDialogWithPlot Creates a modal input dialog with a plot above.
%
% [userValue] = inputDialogWithPlot()
%
% Opens a dialog window containing:
%   - An Axes for plotting.
%   - A Numeric Edit Field for user input.
%   - OK and Cancel buttons.
%
% The user can enter a numeric value, and upon clicking OK, the function
% plots Y = X^2 at the entered X value and returns the X value.
% Clicking Cancel or closing the dialog returns an empty array.

% Initialize output
userValue = [];

% Get screen size
screenSize = get(0, 'ScreenSize'); % [left bottom width height]
dlgWidth = 700;
dlgHeight = 450;

% Calculate position to center the dialog
dlgPosX = (screenSize(3) - dlgWidth) / 2;
dlgPosY = (screenSize(4) - dlgHeight) / 2;

% Create the dialog window as a uifigure
% Make the dialog modal to prevent interaction with other windows
% by using 'modal'
dlg = uifigure('Name', 'Input Threshold Value for Visualization', ...
    'Position', [dlgPosX dlgPosY dlgWidth dlgHeight], ...
    'WindowStyle','modal',...
    'CloseRequestFcn', @(src, event) onCloseRequest(src, event));



% Create axes for plotting within the dialog
ax = uiaxes('Parent', dlg, ...
    'Position', [50, 150, 600, 250], ...
    'Box', 'on');

% Create GMM fit
flattened_data = nonzeros(data(:));

% Compute the histogram normalized to a PDF
[counts, bin_edges] = histcounts(flattened_data, 'Normalization', 'pdf'); % Histogram as PDF
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2; % Bin centers

% Step 2: Define the Lognormal Mixture Model
lognormal_pdf = @(x, mu, sigma) (1 ./ (x * sigma * sqrt(2 * pi))) .* exp(-((log(x) - mu).^2) ./ (2 * sigma^2));
lognormal_mixture = @(params, x) ...
    params(1) * lognormal_pdf(x, params(2), params(3)) + ...
    params(4) * lognormal_pdf(x, params(5), params(6));

% Step 3: Initial Guess for Parameters
% [weight1, mu1, sigma1, weight2, mu2, sigma2]
k = 1.5;
initial_params = [0.5, log(mean(flattened_data)), 0.5, 0.5, log(mean(flattened_data) * k), 0.5];

% Step 4: Fit the Lognormal Mixture Model
% Use lsqcurvefit to optimize the parameters
options = optimoptions('lsqcurvefit', 'Display', 'off');
fitted_params = lsqcurvefit(@(params, x) lognormal_mixture(params, x), ...
                            initial_params, bin_centers, counts, [], [], options);

% Extract the fitted parameters
weight1 = fitted_params(1);
mu1 = fitted_params(2);
sigma1 = fitted_params(3);
weight2 = fitted_params(4);
mu2 = fitted_params(5);
sigma2 = fitted_params(6);

% Step 5: Compute Individual and Combined PDFs
pdf1 = weight1 * lognormal_pdf(bin_centers, mu1, sigma1); % Lognormal Component 1
pdf2 = weight2 * lognormal_pdf(bin_centers, mu2, sigma2); % Lognormal Component 2
pdf_combined = pdf1 + pdf2; % Combined PDF

% Step 6: Find Intersection Points of Lognormal Components
% Define the difference function for intersection
f_diff = @(x) weight1 * lognormal_pdf(x, mu1, sigma1) - weight2 * lognormal_pdf(x, mu2, sigma2);

% Find the intersection points using fzero
x_vals = linspace(min(flattened_data), max(flattened_data), length(flattened_data)); % Fine grid for intersection
intersection_points = [];
for i = 1:length(x_vals) - 1
    if sign(f_diff(x_vals(i))) ~= sign(f_diff(x_vals(i + 1))) % Detect sign change
        intersection_points = [intersection_points, fzero(f_diff, [x_vals(i), x_vals(i + 1)])]; % Refine root
    end
end

% Evaluate PDFs at intersection points
pdf_values_at_intersection = lognormal_mixture(fitted_params, intersection_points);

% Step 7: Plot the Histogram and Fitted PDFs
hold(ax, 'on')
% Plot the histogram
histogram(ax,flattened_data, 'Normalization', 'pdf','NumBins',150);
% Plot individual lognormal components
plot(ax,bin_centers, pdf1, 'r--', 'LineWidth', 2);
plot(ax,bin_centers, pdf2, 'b--', 'LineWidth', 2);
% Plot the combined PDF
plot(ax,bin_centers, pdf_combined, 'g-', 'LineWidth', 2);
% Highlight intersection points
plot(ax,intersection_points, pdf_values_at_intersection, 'go', 'MarkerSize', 8, 'LineWidth', 2);

% Add labels, legend, and title
legend(ax,'Histogram', 'Lognormal Component 1', 'Lognormal Component 2', 'Combined Lognormal Mixture', 'Intersection Points');
title(ax,'Histogram Fit with Two Lognormal Functions and Intersection Points');
xlabel(ax,'Pixel Values (CoV Values)');
ylabel(ax,'Probability Density');
hold(ax, 'off')

% Add a label for the input field
lbl = uilabel(dlg, ...
    'Text', q, ...
    'Position', [80, 100, 450, 22], ...
    'HorizontalAlignment', 'left');

% Add a numeric edit field for user input
editField = uieditfield(dlg, 'numeric', ...
    'Position', [590, 100, 50, 22]);%, ...
% 'ValueChangedFcn', @(src, event) onValueChanged(src, event, ax));

% Add an OK button
btnOK = uibutton(dlg, 'push', ...
    'Text', 'OK', ...
    'Position', [150, 40, 80, 30], ...
    'ButtonPushedFcn', @(src, event) onOK(src, event, dlg, editField, ax));

% Add a Cancel button
btnCancel = uibutton(dlg, 'push', ...
    'Text', 'Cancel', ...
    'Position', [450, 40, 80, 30], ...
    'ButtonPushedFcn', @(src, event) onCancel(src, event, dlg));

% Set KeyPressFcn for the dialog to handle Enter key presses
dlg.KeyPressFcn = @(src, event) onDialogKeyPress(src, event, editField, btnOK, ax);

% Wait for the dialog to close before returning control to the caller
uiwait(dlg);

% Callback for the OK button
    function onOK(src, event, dialogHandle, editHandle, axesHandle)
        % Retrieve the user input
        x = editHandle.Value

        % Validate the input
        if isempty(x) || isnan(x)
            uialert(dialogHandle, 'Please enter a valid numeric value.', 'Invalid Input');
            return;
        end

        % Assign the input value to the output
        userValue = x;

        % Resume execution and close the dialog
        uiresume(dialogHandle);
        delete(dialogHandle);
    end

% Callback for the Cancel button
    function onCancel(src, event, dialogHandle)
        % Assign empty to the output
        userValue = [];

        % Resume execution and close the dialog
        uiresume(dialogHandle);
        delete(dialogHandle);
    end

% Callback for handling dialog close requests (e.g., clicking 'X')
    function onCloseRequest(src, event)
        % Treat as a cancel action
        onCancel([], [], src);
    end

% Callback for key presses within the dialog
    function onDialogKeyPress(src, event, editHandle, okButtonHandle, axesHandle)
        % Check if the pressed key is Enter or Return
        if strcmp(event.Key, 'return') || strcmp(event.Key, 'enter')
            % Check if the edit field is currently focused
            if isequal(src.CurrentObject, editHandle)
                % Trigger the OK button callback
                onOK([], [], src, editHandle, axesHandle);
            else
                onOK([], [], src, editHandle, axesHandle);
                % Optional: Define behavior when Enter is pressed outside the edit field
                % For example, you might want to trigger OK or ignore
            end
        end
    end

end
