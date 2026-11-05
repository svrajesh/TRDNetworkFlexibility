function pairwiseStructureCorrs = structureCorrelations(allh,allJ)
pairwiseStructureCorrs = zeros(size(allh,1),size(allh,1));
for i = 1:size(allh, 1)
    rec1h = allh(i,:);
    rec1J = squeeze(allJ(i,:,:));
    for j = i:size(allh, 1)
        rec2h = allh(j,:);
        rec2J = squeeze(allJ(j,:,:));

        %%
        upperTriRec1J = rec1J(triu(true(size(rec1J)), 1)); % Get upper triangle excluding diagonal
        structure1 = [rec1h, upperTriRec1J']; % Concatenate rec1h with upper triangle of rec1J
        upperTriRec2J = rec2J(triu(true(size(rec2J)), 1)); % Get upper triangle excluding diagonal
        structure2 = [rec2h, upperTriRec2J']; % Concatenate rec1h with upper triangle of rec1J


        pairwiseStructureCorrs(i, j) = corr(structure1(:), structure2(:));
        pairwiseStructureCorrs(j, i) = pairwiseStructureCorrs(i, j);
    end
end






end