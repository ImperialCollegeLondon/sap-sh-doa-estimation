function[Ucell,Scell] = svd_for_cells(R)

Ucell = cell(size(R));
Scell = cell(size(R));
for ii = 1:numel(R)
    if ~isempty(R{ii})
        [Ucell{ii},tmp,~] = svd(R{ii});
        Scell{ii} = diag(tmp);
    end
end