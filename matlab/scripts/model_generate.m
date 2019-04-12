% Generates the boundary value problem for a given mesh
%
%   ratio: ratio in `y`-direction, as cutoff for selecting boundaries
%
function model_generate(inputmesh,outputmat,ratio)
    status = mkdir('../../data/models');
    [V,T,F] = readMESH(inputmesh);
    b = boundary_vertices(V,2,ratio);
    b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
    load = [0; -9.8; 0];
    [~,meshname,~] = fileparts(inputmesh);
    fprintf('Mesh (%s)\n#V: %d\n#T: %d\n#b: %d\n', meshname, size(V,1), size(T,1), size(b,1));
    save(outputmat, ...
        'meshname', ...
        'V', ...
        'T', ...
        'F', ...
        'b', ...
        'load');
end