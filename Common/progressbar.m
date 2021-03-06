classdef progressbar < handle
    properties(Access = protected)
        panel           % Panel on which everything sits
        range_ax        % The progress range axes
        pbar            % The bar representing progress (patch)
        ptext           % Percentage label
        width;
        height;
    end
    properties(Access = public, Dependent = true)
        range           % Progress range
        pvalue          % Current value
        percent         % Percentage complete (relative within range)
        position        % Position of the object (panel)
    end
%     properties(Constant = true)
%         width = 70;
%         height = 20;
%     end
    methods
        % Initializer
        function obj = progressbar(fig, pos, range)
            if nargin < 3
                range = [0 1];
            end
            obj.width = pos(3);
            obj.height = pos(4);
            obj.panel = uipanel('Parent', fig, 'Units', 'pixels', 'Position', pos);
            obj.range_ax = axes('Parent', obj.panel, ...
                'Units', 'pixels', 'Position', [0 0 obj.width obj.height], ...
                'XTickLabel', '', 'XTick', [], 'YTickLabel', '', 'YTick', []);
            obj.pbar = patch([range(1) range(1) range(1) range(1)], [0 0 2 2], ...
                [.75 .75 .9], 'Parent', obj.range_ax);
            obj.ptext = text(obj.width/2 - 20, obj.height/2, '0%', 'Parent', obj.range_ax, ...
                'FontWeight', 'bold', 'Units', 'pixels');
            obj.range = range;
        end

        % Property Access Methods
        function set.range(obj, value)
            % Instead of replotting, just reset the XLim to the
            % extremities of the input range.
            set(obj.range_ax, 'XLim', value([1,end]), 'YLim', [0 2]);
            % Reset progress.
            obj.pvalue = value(1);
        end
        function value = get.range(obj)
            value = get(obj.range_ax, 'XLim');
        end
        function set.pvalue(obj, value)
            % Expects a single value to represent progress value and
            % constructs the selection rectangle from that. If multiple
            % values are passed in, all are ignored but the last, since the
            % left edge of the bar is always the first element of the
            % range.
            set(obj.pbar, 'XData', [obj.range(1) value(end) value(end) obj.range(1)]);
            set(obj.ptext, 'String', sprintf('%3.0f%%', obj.percent * 100));
            drawnow;
        end     
        function value = get.pvalue(obj)
            % The progress bar is actually 2D, but we treat as if it is 1D.
            % Hence the XData is actually an array of four values but we
            % only consider the second (progress maximum).
            limits = get(obj.pbar, 'XData');
            value = limits(2);
        end
        function set.percent(obj, value)
            % Expects a single value between 0 and 1.
            limits = obj.range;
            obj.pvalue = value * (limits(2) - limits(1)) + limits(1);
        end     
        function value = get.percent(obj)            
            limits = obj.range;
            pvalue = obj.pvalue;
            value = (pvalue - limits(1)) / (limits(2) - limits(1));
        end
        function set.position(obj, value)
            % Currently you cannot change the width and height of the
            % object, just the left/bottom position.
            set(obj.panel, 'Position', [value(1) value(2) obj.width obj.height]);
        end
        function value = get.position(obj)
            value = get(obj.panel, 'Position');
        end
    end
end