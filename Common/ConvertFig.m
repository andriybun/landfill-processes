function ConvertFig(srcDir, targetFormat)
    % Function for converting all figures in a given folder from *.fig format to a target format.
    % Parameters:
    %   Source directory - directory containing all *.fig files to be converted
    %   Target format - format to which figures will be converted (default - png)

    if nargin < 2
        targetFormat = '.png';
    end
    
    % Searching using wildcards
    figFiles = dir(fullfile(srcDir, '*.fig'));

    for i = 1:length(figFiles)
        fileName = fullfile(srcDir, figFiles(i).name);
        figH = open(fileName);
%         outerPosition = get(gcf, 'Position');
%         set(gcf, 'Position', outerPosition);
        set(gcf,'color','w');
        frame = getframe(figH);
        [x, ~] = frame2im(frame);
        imwrite(x, strcat(fileName(1:end-4), targetFormat));
        
%         saveas(gcf, strcat(fileName(1:end-4), targetFormat));
        close all;
    end
end