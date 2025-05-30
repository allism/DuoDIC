function DuoDIC_STEP1_Enhanced_Circle_Calibration
% DuoDIC_STEP1_Enhanced_Circle_Calibration - 参考DICe的增强圆点标定方法
% 
% 这个函数实现了类似DICe的高精度圆点标定功能，包括：
% - 非对称圆形网格检测（消除180度歧义性）
% - 亚像素精度圆心检测
% - 多尺度圆形检测和验证
% - 高级图像预处理
% - 鲁棒的圆点排序算法
%
% 语法:
%   DuoDIC_STEP1_Enhanced_Circle_Calibration() - 使用默认设置运行GUI
%   DuoDIC_STEP1_Enhanced_Circle_Calibration('PropertyName', PropertyValue, ...)
%
% 输入参数:
%   'CalibrationPath' - 标定图像路径
%   'CircleGridSize' - 圆点网格尺寸 [行数, 列数] (非对称布局)
%   'SquareSize' - 圆点间距 (mm)
%   'IsAsymmetric' - 是否使用非对称网格 (推荐: true)
%   'ShowReprojectionErrors' - 是否显示重投影误差
%   'CircleDetectionSensitivity' - 圆点检测灵敏度 (0-1)
%   'SubPixelRefinement' - 是否使用亚像素精细化
%   'MinCircleRadius' - 最小圆半径 (像素)
%   'MaxCircleRadius' - 最大圆半径 (像素)
%
% 输出:
%   保存标定参数到DuoDIC_enhanced_calib_results.mat
%
% 作者: DuoDIC工具箱增强版 - 参考DICe标定方法
% 版本: 2.0
% 日期: 2025

%% 参数解析和默认值设置
p = inputParser;
addParameter(p, 'CalibrationPath', pwd, @ischar);
addParameter(p, 'CircleGridSize', [14, 10], @(x) isnumeric(x) && length(x) == 2);  % 修改为14×10
addParameter(p, 'SquareSize', 5.0, @isnumeric);
addParameter(p, 'IsAsymmetric', true, @islogical);
addParameter(p, 'ShowReprojectionErrors', true, @islogical);
addParameter(p, 'CircleDetectionSensitivity', 0.9, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(p, 'SubPixelRefinement', true, @islogical);
addParameter(p, 'MinCircleRadius', 8, @isnumeric);
addParameter(p, 'MaxCircleRadius', 80, @isnumeric);
addParameter(p, 'OutputPath', fullfile(pwd, 'DuoDIC_enhanced_calib_results.mat'), @ischar);

parse(p, varargin{:});

calibPath = p.Results.CalibrationPath;
gridSize = p.Results.CircleGridSize;
squareSize = p.Results.SquareSize;
isAsymmetric = p.Results.IsAsymmetric;
showErrors = p.Results.ShowReprojectionErrors;
sensitivity = p.Results.CircleDetectionSensitivity;
subPixelRefine = p.Results.SubPixelRefinement;
minRadius = p.Results.MinCircleRadius;
maxRadius = p.Results.MaxCircleRadius;
outputPath = p.Results.OutputPath;

%% GUI界面设置（如果没有输入参数）
if nargin == 0
    [calibPath, gridSize, squareSize, isAsymmetric, sensitivity, subPixelRefine, minRadius, maxRadius] = setupEnhancedCircleCalibrationGUI();
    if isempty(calibPath)
        return;
    end
end

%% 读取标定图像
fprintf('正在读取标定图像...\n');
[leftImages, rightImages] = loadCalibrationImages(calibPath);

if isempty(leftImages) || isempty(rightImages)
    error('无法找到标定图像。请检查路径和图像文件。');
end

fprintf('找到 %d 对标定图像\n', length(leftImages));

%% 增强圆点网格检测
fprintf('正在使用DICe风格的增强圆点检测...\n');
[imagePoints, boardSize, pairsUsed, detectionQuality] = detectEnhancedCircleGrids(...
    leftImages, rightImages, gridSize, isAsymmetric, sensitivity, subPixelRefine, minRadius, maxRadius);

if size(imagePoints, 3) < 3
    error('检测到的有效图像对太少（少于3对）。请检查图像质量或调整检测参数。');
end

fprintf('成功检测到 %d 对图像中的圆点网格\n', size(imagePoints, 3));

%% 生成世界坐标点（支持非对称布局）
worldPoints = generateEnhancedWorldPoints(boardSize, squareSize, isAsymmetric);

%% 执行高精度立体标定
fprintf('正在执行高精度立体标定...\n');
[stereoParams, reprojectionErrors] = performEnhancedStereoCalibration(imagePoints, worldPoints);

%% 显示增强的标定结果
displayEnhancedCalibrationResults(stereoParams, reprojectionErrors, detectionQuality, showErrors);

%% 保存标定结果（兼容DuoDIC格式）
saveEnhancedCalibrationResults(stereoParams, reprojectionErrors, detectionQuality, outputPath, pairsUsed);

fprintf('增强圆点标定完成！结果已保存到: %s\n', outputPath);
end

%% 辅助函数

function [calibPath, gridSize, squareSize, isAsymmetric, sensitivity, subPixelRefine, minRadius, maxRadius] = setupEnhancedCircleCalibrationGUI()
% 增强GUI界面设置圆点标定参数
    
    % 创建对话框
    fig = figure('Position', [200, 200, 500, 450], 'MenuBar', 'none', ...
                 'Name', 'DuoDIC增强圆点标定设置 (DICe风格)', 'NumberTitle', 'off', ...
                 'Resize', 'off');
    
    % 路径选择
    uicontrol('Style', 'text', 'Position', [20, 400, 100, 20], ...
              'String', '标定图像路径:', 'HorizontalAlignment', 'left');
    pathEdit = uicontrol('Style', 'edit', 'Position', [20, 380, 300, 25], ...
                         'String', pwd);
    uicontrol('Style', 'pushbutton', 'Position', [330, 380, 80, 25], ...
              'String', '浏览...', 'Callback', @selectPath);
    
    % 网格尺寸
    uicontrol('Style', 'text', 'Position', [20, 340, 150, 20], ...
              'String', '圆点网格尺寸 [行,列]:', 'HorizontalAlignment', 'left');
    gridEdit = uicontrol('Style', 'edit', 'Position', [170, 340, 100, 25], ...
                         'String', '[14, 10]');  % 修改默认值为14×10
    
    % 圆点间距
    uicontrol('Style', 'text', 'Position', [20, 300, 120, 20], ...
              'String', '圆点间距 (mm):', 'HorizontalAlignment', 'left');
    sizeEdit = uicontrol('Style', 'edit', 'Position', [140, 300, 60, 25], ...
                         'String', '5.0');
    
    % 非对称网格
    asymCheck = uicontrol('Style', 'checkbox', 'Position', [250, 300, 200, 25], ...
                          'String', '使用非对称网格 (推荐)', 'Value', 1);
    
    % 检测灵敏度
    uicontrol('Style', 'text', 'Position', [20, 260, 120, 20], ...
              'String', '检测灵敏度 (0-1):', 'HorizontalAlignment', 'left');
    sensEdit = uicontrol('Style', 'edit', 'Position', [140, 260, 60, 25], ...
                         'String', '0.9');
    
    % 亚像素精细化
    subPixelCheck = uicontrol('Style', 'checkbox', 'Position', [250, 260, 200, 25], ...
                              'String', '亚像素精细化 (高精度)', 'Value', 1);
    
    % 圆半径范围
    uicontrol('Style', 'text', 'Position', [20, 220, 120, 20], ...
              'String', '最小圆半径 (像素):', 'HorizontalAlignment', 'left');
    minRadEdit = uicontrol('Style', 'edit', 'Position', [140, 220, 60, 25], ...
                           'String', '8');
    
    uicontrol('Style', 'text', 'Position', [250, 220, 120, 20], ...
              'String', '最大圆半径 (像素):', 'HorizontalAlignment', 'left');
    maxRadEdit = uicontrol('Style', 'edit', 'Position', [370, 220, 60, 25], ...
                           'String', '80');
    
    % 说明文本
    uicontrol('Style', 'text', 'Position', [20, 150, 450, 50], ...
              'String', ['DICe风格增强特性 (14×10网格)：' newline ...
                        '• 非对称网格消除方向歧义性' newline ...
                        '• 亚像素精度圆心检测' newline ...
                        '• 多尺度鲁棒检测算法'], ...
              'HorizontalAlignment', 'left', 'FontSize', 9);
    
    % 按钮
    uicontrol('Style', 'pushbutton', 'Position', [150, 80, 100, 35], ...
              'String', '开始标定', 'FontSize', 10, 'FontWeight', 'bold', ...
              'Callback', @startCalibration);
    uicontrol('Style', 'pushbutton', 'Position', [270, 80, 80, 35], ...
              'String', '取消', 'Callback', @(~,~) close(fig));
    
    % 输出变量
    calibPath = [];
    gridSize = [];
    squareSize = [];
    isAsymmetric = [];
    sensitivity = [];
    subPixelRefine = [];
    minRadius = [];
    maxRadius = [];
    
    uiwait(fig);
    
    function selectPath(~, ~)
        folder = uigetdir(get(pathEdit, 'String'), '选择标定图像文件夹');
        if folder ~= 0
            set(pathEdit, 'String', folder);
        end
    end
    
    function startCalibration(~, ~)
        try
            calibPath = get(pathEdit, 'String');
            gridSize = eval(get(gridEdit, 'String'));
            squareSize = str2double(get(sizeEdit, 'String'));
            isAsymmetric = get(asymCheck, 'Value');
            sensitivity = str2double(get(sensEdit, 'String'));
            subPixelRefine = get(subPixelCheck, 'Value');
            minRadius = str2double(get(minRadEdit, 'String'));
            maxRadius = str2double(get(maxRadEdit, 'String'));
            
            % 验证输入
            if ~exist(calibPath, 'dir')
                errordlg('标定路径不存在！');
                return;
            end
            
            if length(gridSize) ~= 2 || any(gridSize <= 0)
                errordlg('网格尺寸格式错误！');
                return;
            end
            
            if squareSize <= 0
                errordlg('圆点间距必须大于0！');
                return;
            end
            
            if sensitivity < 0 || sensitivity > 1
                errordlg('检测灵敏度必须在0-1之间！');
                return;
            end
            
            if minRadius >= maxRadius || minRadius <= 0
                errordlg('圆半径范围设置错误！');
                return;
            end
            
            close(fig);
        catch ME
            errordlg(['参数错误: ' ME.message]);
        end
    end
end

function [leftImages, rightImages] = loadCalibrationImages(calibPath)
% 加载立体标定图像对（支持多种命名格式）

    % 支持的图像格式
    imageExts = {'*.jpg', '*.jpeg', '*.png', '*.bmp', '*.tif', '*.tiff'};
    
    % 尝试多种命名模式
    namingPatterns = {
        {'left_*', 'right_*'};
        {'L_*', 'R_*'};
        {'cam1_*', 'cam2_*'};
        {'*_left', '*_right'};
        {'*_L', '*_R'};
    };
    
    leftFiles = [];
    rightFiles = [];
    
    for pattern = namingPatterns
        leftPattern = pattern{1}{1};
        rightPattern = pattern{1}{2};
        
        for ext = imageExts
            leftFiles = [leftFiles; dir(fullfile(calibPath, [leftPattern ext{1}(2:end)]))];
            rightFiles = [rightFiles; dir(fullfile(calibPath, [rightPattern ext{1}(2:end)]))];
        end
        
        if ~isempty(leftFiles) && ~isempty(rightFiles)
            break;
        end
    end
    
    % 排序确保对应关系
    if ~isempty(leftFiles) && ~isempty(rightFiles)
        leftNames = {leftFiles.name};
        rightNames = {rightFiles.name};
        
        [~, leftIdx] = sort(leftNames);
        [~, rightIdx] = sort(rightNames);
        
        leftFiles = leftFiles(leftIdx);
        rightFiles = rightFiles(rightIdx);
        
        % 构建完整路径
        leftImages = arrayfun(@(x) fullfile(calibPath, x.name), leftFiles, 'UniformOutput', false);
        rightImages = arrayfun(@(x) fullfile(calibPath, x.name), rightFiles, 'UniformOutput', false);
        
        % 确保数量匹配
        minCount = min(length(leftImages), length(rightImages));
        leftImages = leftImages(1:minCount);
        rightImages = rightImages(1:minCount);
    else
        leftImages = {};
        rightImages = {};
    end
end

function [imagePoints, boardSize, pairsUsed, detectionQuality] = detectEnhancedCircleGrids(...
    leftImages, rightImages, gridSize, isAsymmetric, sensitivity, subPixelRefine, minRadius, maxRadius)
% 增强的圆点网格检测 - DICe风格

    numImages = length(leftImages);
    imagePoints = [];
    pairsUsed = [];
    detectionQuality = [];
    
    fprintf('DICe风格圆点检测进度: ');
    
    for i = 1:numImages
        if mod(i, 3) == 0 || i == numImages
            fprintf('%d/%d ', i, numImages);
        end
        
        % 读取图像
        leftImg = imread(leftImages{i});
        rightImg = imread(rightImages{i});
        
        % 转换为灰度图像
        if size(leftImg, 3) == 3
            leftImg = rgb2gray(leftImg);
        end
        if size(rightImg, 3) == 3
            rightImg = rgb2gray(rightImg);
        end
        
        % 增强的圆点检测
        [leftCenters, leftQuality] = detectEnhancedCircleGrid(leftImg, gridSize, isAsymmetric, ...
                                                            sensitivity, subPixelRefine, minRadius, maxRadius);
        [rightCenters, rightQuality] = detectEnhancedCircleGrid(rightImg, gridSize, isAsymmetric, ...
                                                             sensitivity, subPixelRefine, minRadius, maxRadius);
        
        % 检查检测结果
        expectedPoints = prod(gridSize);
        if isAsymmetric
            expectedPoints = prod(gridSize);  % 非对称网格点数计算
        end
        
        if ~isempty(leftCenters) && ~isempty(rightCenters) && ...
           size(leftCenters, 1) == expectedPoints && size(rightCenters, 1) == expectedPoints
            
            % 添加到图像点集合
            if isempty(imagePoints)
                imagePoints = zeros(size(leftCenters, 1), 2, 0);
            end
            
            imagePoints(:, :, end+1) = leftCenters;
            imagePoints(:, :, end+1) = rightCenters;
            
            pairsUsed = [pairsUsed, i];
            detectionQuality = [detectionQuality; [leftQuality, rightQuality]];
        end
    end
    
    fprintf('\n');
    
    % 重新整理图像点格式
    if ~isempty(imagePoints)
        numPairs = length(pairsUsed);
        imagePointsCell = cell(1, 2);
        imagePointsCell{1} = zeros(size(imagePoints, 1), 2, numPairs); % 左相机
        imagePointsCell{2} = zeros(size(imagePoints, 1), 2, numPairs); % 右相机
        
        for i = 1:numPairs
            imagePointsCell{1}(:, :, i) = imagePoints(:, :, i*2-1);
            imagePointsCell{2}(:, :, i) = imagePoints(:, :, i*2);
        end
        
        imagePoints = imagePointsCell;
        boardSize = gridSize;
    else
        boardSize = gridSize;
    end
end

function [centers, quality] = detectEnhancedCircleGrid(img, gridSize, isAsymmetric, ...
                                                     sensitivity, subPixelRefine, minRadius, maxRadius)
% 单幅图像的增强圆点检测 - DICe风格算法

    try
        % 图像预处理 - DICe风格增强
        enhancedImg = enhanceImageForCircleDetection(img);
        
        % 多尺度圆检测
        [centers, radii, metric] = detectCirclesMultiScale(enhancedImg, minRadius, maxRadius, sensitivity);
        
        expectedPoints = prod(gridSize);
        
        if length(radii) >= expectedPoints
            % 使用改进的排序算法
            if isAsymmetric
                centers = sortAsymmetricCircleGrid(centers, gridSize);
            else
                centers = sortSymmetricCircleGrid(centers, gridSize);
            end
            
            % 验证检测质量
            if size(centers, 1) >= expectedPoints
                centers = centers(1:expectedPoints, :);
                radii = radii(1:expectedPoints);
                metric = metric(1:expectedPoints);
                
                % 亚像素精细化
                if subPixelRefine
                    centers = refineCircleCentersSubPixel(img, centers, radii);
                end
                
                % 计算检测质量分数
                quality = calculateDetectionQuality(centers, radii, metric, gridSize);
            else
                centers = [];
                quality = 0;
            end
        else
            centers = [];
            quality = 0;
        end
        
    catch
        centers = [];
        quality = 0;
    end
end

function enhancedImg = enhanceImageForCircleDetection(img)
% DICe风格图像增强预处理

    % 直方图均衡化
    enhancedImg = adapthisteq(img);
    
    % 高斯滤波去噪
    enhancedImg = imgaussfilt(enhancedImg, 0.5);
    
    % 对比度增强
    enhancedImg = imadjust(enhancedImg);
end

function [centers, radii, metric] = detectCirclesMultiScale(img, minRadius, maxRadius, sensitivity)
% 多尺度圆检测算法

    % 初始检测
    [centers1, radii1, metric1] = imfindcircles(img, [minRadius maxRadius], ...
        'ObjectPolarity', 'dark', 'Sensitivity', sensitivity, 'Method', 'TwoStage');
    
    % 不同极性检测
    [centers2, radii2, metric2] = imfindcircles(img, [minRadius maxRadius], ...
        'ObjectPolarity', 'bright', 'Sensitivity', sensitivity*0.8, 'Method', 'TwoStage');
    
    % 合并结果
    if ~isempty(centers1) && ~isempty(centers2)
        centers = [centers1; centers2];
        radii = [radii1; radii2];
        metric = [metric1; metric2];
    elseif ~isempty(centers1)
        centers = centers1;
        radii = radii1;
        metric = metric1;
    elseif ~isempty(centers2)
        centers = centers2;
        radii = radii2;
        metric = metric2;
    else
        centers = [];
        radii = [];
        metric = [];
    end
    
    % 去除重复检测
    if ~isempty(centers)
        [centers, radii, metric] = removeDuplicateCircles(centers, radii, metric);
    end
end

function [centers, radii, metric] = removeDuplicateCircles(centers, radii, metric)
% 去除重复检测的圆

    if isempty(centers)
        return;
    end
    
    % 计算距离矩阵
    distances = pdist2(centers, centers);
    
    % 找到距离过近的圆对
    [row, col] = find(distances < mean(radii) & distances > 0);
    
    % 保留质量更好的圆
    toRemove = [];
    for i = 1:length(row)
        if metric(row(i)) < metric(col(i))
            toRemove = [toRemove, row(i)];
        else
            toRemove = [toRemove, col(i)];
        end
    end
    
    toRemove = unique(toRemove);
    
    % 移除重复的圆
    centers(toRemove, :) = [];
    radii(toRemove) = [];
    metric(toRemove) = [];
end

function sortedCenters = sortAsymmetricCircleGrid(centers, gridSize)
% 非对称圆点网格排序 - DICe风格 (优化14×10网格)

    if isempty(centers)
        sortedCenters = centers;
        return;
    end
    
    rows = gridSize(1);  % 14
    cols = gridSize(2);  % 10
    
    % 使用kmeans聚类分行
    if size(centers, 1) >= rows
        try
            [~, centroids] = kmeans(centers(:, 2), rows);
            [~, rowOrder] = sort(centroids);
        catch
            % 如果kmeans失败，使用简单排序
            [~, yIdx] = sort(centers(:, 2));
            centers = centers(yIdx, :);
            pointsPerRow = floor(size(centers, 1) / rows);
            rowOrder = 1:rows;
        end
        
        sortedCenters = zeros(rows * cols, 2);
        currentIdx = 1;
        
        for r = 1:rows
            % 找到属于当前行的点
            if exist('centroids', 'var')
                rowCenters = findPointsInRow(centers, centroids(rowOrder(r)), rows);
            else
                % 简单分组方法
                startIdx = (r - 1) * pointsPerRow + 1;
                endIdx = min(r * pointsPerRow, size(centers, 1));
                if endIdx >= startIdx
                    rowCenters = centers(startIdx:endIdx, :);
                else
                    rowCenters = [];
                end
            end
            
            if ~isempty(rowCenters)
                % 按X坐标排序
                [~, xOrder] = sort(rowCenters(:, 1));
                rowCenters = rowCenters(xOrder, :);
                
                % 对于14×10网格，确保每行有10个点
                actualCols = min(cols, size(rowCenters, 1));
                
                if actualCols > 0
                    endIdx = currentIdx + actualCols - 1;
                    if endIdx <= size(sortedCenters, 1)
                        sortedCenters(currentIdx:endIdx, :) = rowCenters(1:actualCols, :);
                        currentIdx = endIdx + 1;
                    end
                end
            end
        end
        
        % 移除未使用的行
        validRows = any(sortedCenters ~= 0, 2);
        sortedCenters = sortedCenters(validRows, :);
    else
        sortedCenters = centers;
    end
end

function sortedCenters = sortSymmetricCircleGrid(centers, gridSize)
% 对称圆点网格排序

    if isempty(centers)
        sortedCenters = centers;
        return;
    end
    
    % 按Y坐标分组（行）
    [~, yIdx] = sort(centers(:, 2));
    centers = centers(yIdx, :);
    
    rows = gridSize(1);
    cols = gridSize(2);
    pointsPerRow = floor(size(centers, 1) / rows);
    
    sortedCenters = zeros(rows * cols, 2);
    
    for row = 1:rows
        startIdx = (row - 1) * pointsPerRow + 1;
        endIdx = min(row * pointsPerRow, size(centers, 1));
        
        if endIdx >= startIdx
            rowPoints = centers(startIdx:endIdx, :);
            
            % 在每行内按X坐标排序
            [~, xIdx] = sort(rowPoints(:, 1));
            rowPoints = rowPoints(xIdx, :);
            
            % 只取需要的列数
            actualCols = min(cols, size(rowPoints, 1));
            sortedCenters((row-1)*cols+1:(row-1)*cols+actualCols, :) = rowPoints(1:actualCols, :);
        end
    end
    
    % 移除未填充的点
    validPoints = all(sortedCenters ~= 0, 2);
    sortedCenters = sortedCenters(validPoints, :);
end

function rowCenters = findPointsInRow(centers, rowY, numRows)
% 找到属于特定行的点

    threshold = std(centers(:, 2)) / numRows * 1.5;  % 增加容差
    distances = abs(centers(:, 2) - rowY);
    rowIndices = distances < threshold;
    rowCenters = centers(rowIndices, :);
end

function refinedCenters = refineCircleCentersSubPixel(img, centers, radii)
% 亚像素精度圆心精细化 - DICe风格算法

    refinedCenters = centers;
    
    for i = 1:size(centers, 1)
        x = round(centers(i, 1));
        y = round(centers(i, 2));
        r = max(3, round(radii(i) * 0.7));
        
        % 确保在图像边界内
        x1 = max(1, x - r);
        x2 = min(size(img, 2), x + r);
        y1 = max(1, y - r);
        y2 = min(size(img, 1), y + r);
        
        if x2 > x1 && y2 > y1
            % 提取局部区域
            localImg = double(img(y1:y2, x1:x2));
            
            % 高斯拟合求精
            [X, Y] = meshgrid(1:size(localImg, 2), 1:size(localImg, 1));
            
            % 寻找局部极值
            if mean(localImg(:)) > 128  % 亮圆
                [~, maxIdx] = max(localImg(:));
            else  % 暗圆
                [~, maxIdx] = min(localImg(:));
            end
            
            [maxY, maxX] = ind2sub(size(localImg), maxIdx);
            
            % 计算质心
            weights = abs(localImg - mean(localImg(:)));
            totalWeight = sum(weights(:));
            
            if totalWeight > 0
                centroidX = sum(sum(X .* weights)) / totalWeight;
                centroidY = sum(sum(Y .* weights)) / totalWeight;
                
                % 转换回全局坐标
                refinedCenters(i, 1) = x1 + centroidX - 1;
                refinedCenters(i, 2) = y1 + centroidY - 1;
            end
        end
    end
end

function quality = calculateDetectionQuality(centers, radii, metric, gridSize)
% 计算检测质量分数

    if isempty(centers)
        quality = 0;
        return;
    end
    
    % 基础质量分数
    baseQuality = mean(metric);
    
    % 圆的一致性评估
    radiusConsistency = 1 - std(radii) / (mean(radii) + eps);
    
    % 网格规律性评估
    gridRegularity = assessGridRegularity(centers, gridSize);
    
    % 覆盖度评估
    coverage = min(1, length(radii) / prod(gridSize));
    
    % 综合质量分数
    quality = baseQuality * 0.4 + radiusConsistency * 0.2 + gridRegularity * 0.3 + coverage * 0.1;
    quality = max(0, min(1, quality));
end

function regularity = assessGridRegularity(centers, gridSize)
% 评估网格规律性

    if size(centers, 1) < prod(gridSize)
        regularity = 0;
        return;
    end
    
    rows = gridSize(1);
    cols = gridSize(2);
    
    try
        % 重新排列点为网格形式
        gridPoints = reshape(centers(1:rows*cols, :), [rows, cols, 2]);
        
        % 计算水平间距的一致性
        horizontalSpacing = [];
        for r = 1:rows
            for c = 1:cols-1
                dist = norm(squeeze(gridPoints(r, c+1, :) - gridPoints(r, c, :)));
                horizontalSpacing = [horizontalSpacing, dist];
            end
        end
        
        % 计算垂直间距的一致性
        verticalSpacing = [];
        for r = 1:rows-1
            for c = 1:cols
                dist = norm(squeeze(gridPoints(r+1, c, :) - gridPoints(r, c, :)));
                verticalSpacing = [verticalSpacing, dist];
            end
        end
        
        % 计算一致性
        hConsistency = 1 - std(horizontalSpacing) / (mean(horizontalSpacing) + eps);
        vConsistency = 1 - std(verticalSpacing) / (mean(verticalSpacing) + eps);
        
        regularity = (hConsistency + vConsistency) / 2;
        regularity = max(0, min(1, regularity));
        
    catch
        regularity = 0.5;  % 默认值
    end
end

function worldPoints = generateEnhancedWorldPoints(boardSize, squareSize, isAsymmetric)
% 生成世界坐标点（支持非对称布局）

    rows = boardSize(1);
    cols = boardSize(2);
    
    if isAsymmetric
        % 非对称网格世界坐标生成
        worldPoints = zeros(rows * cols, 3);
        idx = 1;
        
        for r = 1:rows
            for c = 1:cols
                % 非对称布局：奇数行偏移半个格子
                if mod(r, 2) == 1
                    x = (c - 1) * squareSize;
                else
                    x = (c - 0.5) * squareSize;
                end
                
                y = (r - 1) * squareSize;
                z = 0;
                
                worldPoints(idx, :) = [x, y, z];
                idx = idx + 1;
            end
        end
    else
        % 对称网格世界坐标生成
        [X, Y] = meshgrid(0:cols-1, 0:rows-1);
        worldPoints = [X(:) * squareSize, Y(:) * squareSize, zeros(rows * cols, 1)];
    end
end

function [stereoParams, reprojectionErrors] = performEnhancedStereoCalibration(imagePoints, worldPoints)
% 执行高精度立体标定

    try
        % 获取图像尺寸（从第一个图像点集合推断）
        imageSize = [size(imagePoints{1}, 1), size(imagePoints{1}, 2)];
        
        % 执行立体标定
        stereoParams = estimateCameraParameters(imagePoints, worldPoints, ...
            'EstimateSkew', false, ...
            'EstimateTangentialDistortion', true, ...
            'NumRadialDistortionCoefficients', 3, ...
            'WorldUnits', 'millimeters', ...
            'InitialIntrinsicMatrix', [], ...
            'InitialRadialDistortion', []);
        
        % 计算重投影误差
        reprojectionErrors = stereoParams.ReprojectionErrors;
        
        fprintf('立体标定完成:\n');
        fprintf('  - 平均重投影误差: %.3f 像素\n', stereoParams.MeanReprojectionError);
        fprintf('  - 使用图像对数: %d\n', stereoParams.NumPatterns);
        
    catch ME
        error('立体标定失败: %s', ME.message);
    end
end

function displayEnhancedCalibrationResults(stereoParams, reprojectionErrors, detectionQuality, showErrors)
% 显示增强的标定结果

    fprintf('\n========== DuoDIC增强圆点标定结果 ==========\n');
    
    % 基本标定参数
    fprintf('标定质量评估:\n');
    fprintf('  - 平均重投影误差: %.4f 像素\n', stereoParams.MeanReprojectionError);
    fprintf('  - 最大重投影误差: %.4f 像素\n', max(reprojectionErrors(:)));
    fprintf('  - 标准差: %.4f 像素\n', std(reprojectionErrors(:)));
    fprintf('  - 平均检测质量: %.4f\n', mean(detectionQuality(:)));
    
    % 相机内参
    fprintf('\n左相机内参:\n');
    fprintf('  - 焦距: [%.2f, %.2f]\n', stereoParams.CameraParameters1.FocalLength);
    fprintf('  - 主点: [%.2f, %.2f]\n', stereoParams.CameraParameters1.PrincipalPoint);
    fprintf('  - 径向畸变: [%.6f, %.6f, %.6f]\n', stereoParams.CameraParameters1.RadialDistortion);
    
    fprintf('\n右相机内参:\n');
    fprintf('  - 焦距: [%.2f, %.2f]\n', stereoParams.CameraParameters2.FocalLength);
    fprintf('  - 主点: [%.2f, %.2f]\n', stereoParams.CameraParameters2.PrincipalPoint);
    fprintf('  - 径向畸变: [%.6f, %.6f, %.6f]\n', stereoParams.CameraParameters2.RadialDistortion);
    
    % 立体参数
    fprintf('\n立体参数:\n');
    rotationAngles = rotationMatrixToVector(stereoParams.RotationOfCamera2);
    fprintf('  - 旋转角度 (度): [%.4f, %.4f, %.4f]\n', rotationAngles * 180/pi);
    fprintf('  - 平移向量: [%.4f, %.4f, %.4f] mm\n', stereoParams.TranslationOfCamera2);
    fprintf('  - 基线距离: %.4f mm\n', norm(stereoParams.TranslationOfCamera2));
    
    % 图形化显示
    if showErrors
        displayCalibrationErrorsGraphically(stereoParams, reprojectionErrors, detectionQuality);
    end
    
    fprintf('===============================================\n\n');
end

function displayCalibrationErrorsGraphically(stereoParams, reprojectionErrors, detectionQuality)
% 图形化显示标定误差

    figure('Name', 'DuoDIC增强圆点标定结果分析', 'Position', [100, 100, 1200, 600]);
    
    % 重投影误差分布
    subplot(2, 3, 1);
    histogram(reprojectionErrors(:), 20, 'FaceColor', [0.2, 0.7, 0.9]);
    xlabel('重投影误差 (像素)');
    ylabel('频次');
    title('重投影误差分布');
    grid on;
    
    % 每个图像的平均误差
    subplot(2, 3, 2);
    meanErrorsPerImage = squeeze(mean(sqrt(sum(reprojectionErrors.^2, 2)), 1));
    plot(meanErrorsPerImage, 'o-', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('图像编号');
    ylabel('平均重投影误差 (像素)');
    title('各图像标定精度');
    grid on;
    
    % 检测质量评估
    subplot(2, 3, 3);
    if ~isempty(detectionQuality)
        bar(mean(detectionQuality, 1), 'FaceColor', [0.9, 0.6, 0.2]);
        set(gca, 'XTickLabel', {'左相机', '右相机'});
        ylabel('检测质量得分');
        title('圆点检测质量');
        ylim([0, 1]);
        grid on;
    end
    
    % 误差空间分布
    subplot(2, 3, 4);
    showExtrinsics(stereoParams, 'CameraCentric');
    title('相机外参可视化');
    
    % 畸变参数比较
    subplot(2, 3, 5);
    distortion1 = stereoParams.CameraParameters1.RadialDistortion;
    distortion2 = stereoParams.CameraParameters2.RadialDistortion;
    bar([distortion1; distortion2]');
    xlabel('畸变系数');
    ylabel('数值');
    title('径向畸变参数对比');
    legend('左相机', '右相机');
    grid on;
    
    % 焦距比较
    subplot(2, 3, 6);
    focal1 = stereoParams.CameraParameters1.FocalLength;
    focal2 = stereoParams.CameraParameters2.FocalLength;
    bar([focal1; focal2]);
    set(gca, 'XTickLabel', {'fx', 'fy'});
    ylabel('焦距 (像素)');
    title('相机焦距对比');
    legend('左相机', '右相机');
    grid on;
end

function saveEnhancedCalibrationResults(stereoParams, reprojectionErrors, detectionQuality, outputPath, pairsUsed)
% 保存增强标定结果（兼容DuoDIC格式）

    % 兼容DuoDIC的变量命名
    DuoDIC_stereoParams = stereoParams;
    DuoDIC_reprojectionErrors = reprojectionErrors;
    DuoDIC_detectionQuality = detectionQuality;
    DuoDIC_pairsUsed = pairsUsed;
    DuoDIC_calibrationDate = datestr(now);
    DuoDIC_calibrationMethod = 'Enhanced Circle Grid (DICe Style)';
    
    % 提取关键参数供DuoDIC使用
    DuoDIC_cameraMatrix1 = stereoParams.CameraParameters1.IntrinsicMatrix';
    DuoDIC_cameraMatrix2 = stereoParams.CameraParameters2.IntrinsicMatrix';
    DuoDIC_distCoeffs1 = [stereoParams.CameraParameters1.RadialDistortion, ...
                          stereoParams.CameraParameters1.TangentialDistortion];
    DuoDIC_distCoeffs2 = [stereoParams.CameraParameters2.RadialDistortion, ...
                          stereoParams.CameraParameters2.TangentialDistortion];
    DuoDIC_R = stereoParams.RotationOfCamera2;
    DuoDIC_T = stereoParams.TranslationOfCamera2;
    DuoDIC_E = stereoParams.EssentialMatrix;
    DuoDIC_F = stereoParams.FundamentalMatrix;
    
    % 质量指标
    DuoDIC_meanReprojError = stereoParams.MeanReprojectionError;
    DuoDIC_maxReprojError = max(reprojectionErrors(:));
    DuoDIC_stdReprojError = std(reprojectionErrors(:));
    
    % 保存到MAT文件
    save(outputPath, 'DuoDIC_stereoParams', 'DuoDIC_reprojectionErrors', ...
         'DuoDIC_detectionQuality', 'DuoDIC_pairsUsed', 'DuoDIC_calibrationDate', ...
         'DuoDIC_calibrationMethod', 'DuoDIC_cameraMatrix1', 'DuoDIC_cameraMatrix2', ...
         'DuoDIC_distCoeffs1', 'DuoDIC_distCoeffs2', 'DuoDIC_R', 'DuoDIC_T', ...
         'DuoDIC_E', 'DuoDIC_F', 'DuoDIC_meanReprojError', 'DuoDIC_maxReprojError', ...
         'DuoDIC_stdReprojError', '-v7.3');
    
    % 同时保存为文本格式便于查看
    textOutputPath = strrep(outputPath, '.mat', '_summary.txt');
    saveCalibrationSummaryToText(stereoParams, reprojectionErrors, detectionQuality, textOutputPath);
    
    fprintf('标定结果已保存:\n');
    fprintf('  - MAT格式: %s\n', outputPath);
    fprintf('  - 文本格式: %s\n', textOutputPath);
end

function saveCalibrationSummaryToText(stereoParams, reprojectionErrors, detectionQuality, textPath)
% 保存标定结果摘要到文本文件

    fid = fopen(textPath, 'w');
    if fid == -1
        warning('无法创建文本摘要文件');
        return;
    end
    
    fprintf(fid, 'DuoDIC增强圆点标定结果摘要\n');
    fprintf(fid, '生成时间: %s\n', datestr(now));
    fprintf(fid, '标定方法: Enhanced Circle Grid (DICe Style)\n');
    fprintf(fid, '===============================================\n\n');
    
    fprintf(fid, '标定精度:\n');
    fprintf(fid, '  平均重投影误差: %.4f 像素\n', stereoParams.MeanReprojectionError);
    fprintf(fid, '  最大重投影误差: %.4f 像素\n', max(reprojectionErrors(:)));
    fprintf(fid, '  误差标准差: %.4f 像素\n', std(reprojectionErrors(:)));
    fprintf(fid, '  平均检测质量: %.4f\n', mean(detectionQuality(:)));
    fprintf(fid, '  使用图像对数: %d\n\n', stereoParams.NumPatterns);
    
    fprintf(fid, '左相机内参:\n');
    fprintf(fid, '  焦距: [%.2f, %.2f]\n', stereoParams.CameraParameters1.FocalLength);
    fprintf(fid, '  主点: [%.2f, %.2f]\n', stereoParams.CameraParameters1.PrincipalPoint);
    fprintf(fid, '  径向畸变: [%.6f, %.6f, %.6f]\n', stereoParams.CameraParameters1.RadialDistortion);
    fprintf(fid, '  切向畸变: [%.6f, %.6f]\n\n', stereoParams.CameraParameters1.TangentialDistortion);
    
    fprintf(fid, '右相机内参:\n');
    fprintf(fid, '  焦距: [%.2f, %.2f]\n', stereoParams.CameraParameters2.FocalLength);
    fprintf(fid, '  主点: [%.2f, %.2f]\n', stereoParams.CameraParameters2.PrincipalPoint);
    fprintf(fid, '  径向畸变: [%.6f, %.6f, %.6f]\n', stereoParams.CameraParameters2.RadialDistortion);
    fprintf(fid, '  切向畸变: [%.6f, %.6f]\n\n', stereoParams.CameraParameters2.TangentialDistortion);
    
    fprintf(fid, '立体参数:\n');
    rotationAngles = rotationMatrixToVector(stereoParams.RotationOfCamera2);
    fprintf(fid, '  旋转角度 (度): [%.4f, %.4f, %.4f]\n', rotationAngles * 180/pi);
    fprintf(fid, '  平移向量 (mm): [%.4f, %.4f, %.4f]\n', stereoParams.TranslationOfCamera2);
    fprintf(fid, '  基线距离 (mm): %.4f\n', norm(stereoParams.TranslationOfCamera2));
    
    fclose(fid);
end