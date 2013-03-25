function Merged = MergeStructs(Struct1, Struct2, varargin)
    MergedCell = cat(2, fieldnames(Struct1), struct2cell(Struct1));
    MergedCell = cat(1, MergedCell, cat(2, fieldnames(Struct2), struct2cell(Struct2)));
    
    %% TODO: treat varargin
    %%
    
    if (numel(unique(MergedCell(:, 1))) ~= numel(fieldnames(Struct1)) + numel(fieldnames(Struct2)))
        error('StructMerge:Error', 'Structures have common fields');
    end
    
    Merged = cell2struct(MergedCell(:, 2), MergedCell(:, 1));
end