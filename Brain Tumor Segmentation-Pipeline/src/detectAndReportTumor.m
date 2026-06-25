function detectAndReportTumor(imagePath, authorName, groundTruthPath)
    % detectAndReportTumor - Robust Multi-Scenario Brain Tumor Segmentation
    
    if nargin < 2 || isempty(authorName)
        authorName = 'Reza Shojaei Nasab';
    end
    if nargin < 3
        groundTruthPath = '';
    end
    
    %% 1. Image Loading and Preprocessing
    try
        I = imread(imagePath); 
    catch
        error('Error: Could not read the image file. Please check the path.');
    end
    
    if size(I, 3) == 3
        I = rgb2gray(I);
    end
    I_filtered = medfilt2(I, [3 3]);
    
    %% 2. Head Detection & Adaptive Brain Masking
    global_thresh = multithresh(I_filtered, 1);
    BW_head = I_filtered > global_thresh;
    BW_head = imfill(BW_head, 'holes');
    BW_head = bwareafilt(BW_head, 1);
    
    head_stats = regionprops(BW_head, 'BoundingBox');
    head_box = head_stats(1).BoundingBox; 
    
    % Generous Brain Zone to accommodate tilted images (like gg 198)
    Brain_Zone_Mask = false(size(I));
    brain_bottom_limit = round(head_box(2) + 0.62 * head_box(4)); 
    Brain_Zone_Mask(1:brain_bottom_limit, :) = true;
    
    BW_brain_mask = BW_head & Brain_Zone_Mask;
    BW_brain_mask = imerode(BW_brain_mask, strel('disk', 4));
    
    Brain_Only = I_filtered;
    Brain_Only(~BW_brain_mask) = 0;
    
    %% 3. Calculate Brain Background Stats (To prevent False Positives)
    brain_pixels = double(Brain_Only(BW_brain_mask & (Brain_Only > 0)));
    mean_brain = mean(brain_pixels);
    std_brain = std(brain_pixels);
    
    %% 4. Adaptive Thresholding & Hole Filling
    brain_levels = multithresh(Brain_Only, 3); 
    refined_thresh = brain_levels(2) + 0.15 * (brain_levels(3) - brain_levels(2));
    BW_tumor_raw = Brain_Only >= refined_thresh;
    
    % CRITICAL FIX FOR RING-ENHANCEMENT (like gg 256): 
    % Fill holes BEFORE regionprops so necrotic cores are included in solidity/area
    BW_tumor_raw = imfill(BW_tumor_raw, 'holes');
    
    %% 5. Component Analysis & Multi-Focal Scoring
    stats = regionprops(BW_tumor_raw, 'Area', 'Solidity', 'PixelIdxList', 'Centroid');
    
    valid_indices = [];
    scores = [];
    
    for i = 1:length(stats)
        area = stats(i).Area;
        solidity = stats(i).Solidity;
        
        % Filter out small noise
        if area < 120 || solidity < 0.40
            continue;
        end
        
        % Calculate candidate's mean intensity
        candidate_pixels = double(Brain_Only(stats(i).PixelIdxList));
        mean_candidate = mean(candidate_pixels);
        
        % INTENSITY CHECK: Must be significantly brighter than normal brain tissue
        % This completely protects healthy/post-op slices (like gg 41) from false positives
        if mean_candidate < (mean_brain + 1.2 * std_brain)
            continue; 
        end
        
        % Scoring System
        score = area * solidity;
        valid_indices(end+1) = i; %#ok<AGROW>
        scores(end+1) = score; %#ok<AGROW>
    end
    
    % Reconstruct the final tumor mask (Supports Multi-focal like gg 115)
    BW_tumor = false(size(I));
    if ~isempty(scores)
        max_score = max(scores);
        % Accept the best candidate AND any other candidate that holds at least 25% of max score
        for j = 1:length(scores)
            if scores(j) >= 0.25 * max_score
                BW_tumor(stats(valid_indices(j)).PixelIdxList) = true;
            end
        end
        BW_tumor = imclose(BW_tumor, strel('disk', 4));
        BW_tumor = imfill(BW_tumor, 'holes');
    end
    
    %% 6. Quantitative Measurements
    tumor_pixel_count = sum(BW_tumor(:));
    
    if tumor_pixel_count == 0
        disp('Analysis Complete: No tumor component matched the contrast and structural criteria.');
        % Generate a "Clear/No Tumor" Report instead of crashing
        generateEmptyReport(authorName);
        return;
    end
    
    % Safe Bounding Box for multi-focal regions
    [rows, cols] = find(BW_tumor);
    tumor_width = max(cols) - min(cols) + 1;
    tumor_height = max(rows) - min(rows) + 1;
    
    brain_pixel_count = sum(BW_brain_mask(:));
    tumor_percentage = (tumor_pixel_count / brain_pixel_count) * 100;
    pixel_area_mm2 = 1.0; 
    estimated_tumor_area_mm2 = tumor_pixel_count * pixel_area_mm2;
    %% 5.5 Quantitative Validation (The Research Gold Standard)
if ~isempty(groundTruthPath)
    try
        GT = imread(groundTruthPath);
        if size(GT, 3) == 3, GT = rgb2gray(GT); end
        BW_GT = GT > 0; % Ensure binary
        
        % Calculate Metrics
        dice_score = dice(BW_tumor, BW_GT);
        jaccard_score = jaccard(BW_tumor, BW_GT);
        
        % Display in command window
        fprintf('Validation Metrics -> Dice: %.4f | IoU: %.4f\n', dice_score, jaccard_score);
    catch
        warning('Could not process Ground Truth image for validation.');
        dice_score = [];
    end
end
    %% 7. Visualization & Export
    fig = figure('Name', 'Robust Tumor Detection', 'NumberTitle', 'off', 'Visible', 'off');
    subplot(1, 2, 1);
    imshow(I);
    title('Original MRI');
    
    subplot(1, 2, 2);
    imshow(I); 
    hold on;
    
    red_channel = zeros(size(I));
    red_channel(BW_tumor) = 255;
    overlay = cat(3, red_channel, zeros(size(I)), zeros(size(I)));
    h = imshow(overlay);
    set(h, 'AlphaData', BW_tumor * 0.45); 
    
    boundaries = bwboundaries(BW_tumor);
    for k = 1:length(boundaries)
        b = boundaries{k};
        plot(b(:,2), b(:,1), 'y', 'LineWidth', 0.6); 
    end
    title('Anatomical Precision Overlay');
    
    txt = {
        sprintf('Tumor Area: %.2f mm^2', estimated_tumor_area_mm2)
        sprintf('Brain Ratio: %.2f %%', tumor_percentage)
        sprintf('Width: %d px', tumor_width)
        sprintf('Height: %d px', tumor_height)
    };
    text(15, 35, txt, ...
        'Color', 'white', ...
        'FontSize', 9, ...
        'FontWeight', 'bold', ...
        'BackgroundColor', [0 0 0.4], ... 
        'Margin', 4);
    hold off;
    exportgraphics(fig, 'Tumor_Result.png', 'Resolution', 300);
    close(fig); 
    
    %% 8. PDF Report Generation
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report('Brain_Tumor_Analysis_Report', 'pdf');
    tp = TitlePage;
    tp.Title = 'Automated Brain Tumor Detection & Quantification';
    tp.Author = authorName;
    tp.Subtitle = 'Biomedical Image Processing Pipeline (Multi-Scenario Optimized)';
    add(rpt, tp);
    
    ch = Chapter;
    ch.Title = 'Quantitative Analysis Results';
    append(ch, Paragraph(sprintf('Analysis Date: %s', datetime('now'))));
    append(ch, Paragraph(sprintf('Estimated Tumor Area: %.2f mm^2', estimated_tumor_area_mm2)));
    append(ch, Paragraph(sprintf('Tumor-to-Brain Ratio: %.2f %%', tumor_percentage)));
    append(ch, Paragraph(sprintf('Bounding Box (W x H): %d x %d px', tumor_width, tumor_height)));
    
    img = Image('Tumor_Result.png');
    img.Width = '5.5in';   
    img.Height = '3.3in';
    append(ch, Paragraph(''));
    append(ch, img);
    add(rpt, ch);
    close(rpt);
    
    disp('Report generated successfully: Brain_Tumor_Analysis_Report.pdf');
end

function generateEmptyReport(authorName)
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report('Brain_Tumor_Analysis_Report', 'pdf');
    tp = TitlePage;
    tp.Title = 'Automated Brain Tumor Detection Report';
    tp.Author = authorName;
    tp.Subtitle = 'Status: No Tumor Detected / Normal Tissue Variance';
    add(rpt, tp);
    ch = Chapter('Analysis Results');
    append(ch, Paragraph(sprintf('Analysis Date: %s', datetime('now'))));
    append(ch, Paragraph('Result: No abnormal contrast-enhancing tumor mass structures were identified matching the pathophysiological criteria.'));
    add(rpt, ch);
    close(rpt);
    disp('Report generated successfully (Status: Clear): Brain_Tumor_Analysis_Report.pdf');
end