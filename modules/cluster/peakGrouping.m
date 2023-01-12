function [objOut] = peakGrouping(detMat)

numDetectedObjects = size(detMat,2);
objOut = [];

% sort the detMat matrix according to the cell power
[~, order] = sort(detMat(3,:), 'descend');
detMat = detMat(:,order);

kernel_size_rng = 3;
kernel_size_dop = 3;

for ni = 1:numDetectedObjects
    detectedObjFlag = 1;
    rangeIdx = detMat(2,ni);
    dopplerIdx = detMat(1,ni);
    peakVal = detMat(3,ni);
    kernel = zeros(kernel_size_rng,kernel_size_dop);
    
    %% fill the middle column of the  kernel
    kernel(ceil(kernel_size_rng/2),ceil(kernel_size_dop/2)) = peakVal;
    
    for i=1:kernel_size_rng
        for j=1:kernel_size_dop
            need_index = find(detMat(1,:) == dopplerIdx+j-ceil(kernel_size_dop/2) & detMat(2,:) == rangeIdx+i-ceil(kernel_size_rng/2));
            if ~isempty(need_index)
                kernel(i,j) = detMat(3,need_index(1));
            end
        end
    end
    
%     need_index = find(detMat(1,:) == dopplerIdx & detMat(2,:) == rangeIdx+1);
%     if ~isempty(need_index)
%         kernel(1,2) = detMat(3,need_index(1));
%     end
%     
%     need_index = find(detMat(1,:) == dopplerIdx & detMat(2,:) == rangeIdx-1);
%     if ~isempty(need_index)
%         kernel(3,2) = detMat(3,need_index(1));
%     end
% 
%     % fill the left column of the kernal
%     need_index = find(detMat(1,:) == dopplerIdx-1 & detMat(2,:) == rangeIdx+1);
%     if ~isempty(need_index)
%         kernel(1,1) = detMat(3,need_index(1));
%     end
%     
%     need_index = find(detMat(1,:) == dopplerIdx-1 & detMat(2,:) == rangeIdx);
%     if ~isempty(need_index)
%         kernel(2,1) = detMat(3,need_index(1));
%     end
%     
%     need_index = find(detMat(1,:) == dopplerIdx-1 & detMat(2,:) == rangeIdx-1);
%     if ~isempty(need_index)
%         kernel(3,1) = detMat(3,need_index(1));
%     end
%     
%     % Fill the right column of the kernel
%     need_index = find(detMat(1,:) == dopplerIdx+1 & detMat(2,:) == rangeIdx+1);
%     if ~isempty(need_index)
%         kernel(1,3) = detMat(3,need_index(1));
%     end
%     
%     need_index = find(detMat(1,:) == dopplerIdx+1 & detMat(2,:) == rangeIdx);
%     if ~isempty(need_index)
%         kernel(2,3) = detMat(3,need_index(1));
%     end
%     
%     need_index = find(detMat(1,:) == dopplerIdx+1 & detMat(2,:) == rangeIdx-1);
%     if ~isempty(need_index)
%         kernel(3,3) = detMat(3,need_index(1));
%     end
    
    % Compare the detected object to its neighbors.Detected object is
    % at index [2,2]
    if peakVal ~= max(max(kernel))
        detectedObjFlag = 0;
    end
    
    if detectedObjFlag == 1
        objOut = [objOut, detMat(:,ni)];
    end
end
end