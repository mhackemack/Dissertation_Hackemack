function loadmydata(dataID,porder)

global dat npar snq mats

sn=snq.sn;
mu=snq.mu;

switch dataID
    case 1 % problem A from paper
        % number of elements per zone
        nel_zone = [ 40 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qv=[0.01];
        % incoming flux values
        inc(1:sn)   = 0;

    case 2 % problem B from paper
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qv=[0];
        % incoming flux values
        inc(1:sn)   = 0;
        inc(sn/2+1:end) = 1;

    case 3 % problem C from paper
        % number of elements per zone
        nel_zone = [ 50 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qv=[0];
        % incoming flux values
        inc(1:sn)   = 0;
        inc(sn)=1/SNQ.w(sn);

    case 4 % problem D from paper
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 1 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[1/epsilon];
        sigs=[1/epsilon-epsilon];
        % volumetric source value, per zone
        qv=[epsilon];
        % incoming flux values
        inc(1:sn)   = 0;

    case 5 % 2-zone problem
        % number of elements per zone
        nel_zone = [ 50 50 ];
        % width of each zone
        width_zone = [ 0.5 1 ];
        % sigt/sigs per zone
        sigt=[100 1e4];
        sigs=sigt;
        % volumetric source value, per zone
        qv=[1e-2 1e-2];
        % incoming flux values
        inc(1:sn)   = 0;

    case 6 % 3-zone problem
        % number of elements per zone
        nel_zone = [ 50 50 10]*3;
        % width of each zone
        width_zone = [ 1/3 2/3 1 ];
        % sigt/sigs per zone
        sigt=[100 100 100];
        sigs=[50 50 50]*2;
        % volumetric source value, per zone
        qv=[1 2 1];
        % incoming flux values
        inc(1:sn)   = 0;

    case 7 % Reed (5-zone problem)
        % number of elements per zone
        nel_zone = [ 20 10 20 10 20]*2;
        % width of each zone
        width_zone = [ 2 3 5 6 8 ];
        % sigt/sigs per zone
        sigt=[50 5 1e-8 1 1 ];
        sigs=[ 0 0    0 0.9 0.9];
        % volumetric source value, per zone
        qv=[50 0 0 1 0];
        % incoming flux values
        inc(1:sn)   = 0;

    case 8 % Adams (PIM)
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 500 ];
        % sigt/sigs per zone
        sigt=[1];
        sigs=[1];
        % volumetric source value, per zone
        qv=[0];
        % incoming flux values
        inc(1:sn)   = 0;
%         inc(sn/2+1) = 1;
        inc(sn/2+1)=1/snq.w(sn/2+1);

    case 9 % Adams (PIM) ref
        % number of elements per zone
        nel_zone = [ 1000 1000];
        % width of each zone
        width_zone = [ 5 500 ];
        % sigt/sigs per zone
        sigt=[1 1];
        sigs=[1 1];
        % volumetric source value, per zone
        qv=[0 0];
        % incoming flux values
        inc(1:sn)   = 0;
%         inc(sn/2+1) = 1;
        inc(sn/2+1)=1/snq.w(sn/2+1);
    % fine mesh limit - DSA testing
    case 10
        % number of elements per zone
        nel_zone = [ 1000 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        sigt=[10];
        sigs=[9.999];
        % volumetric source value, per zone
        qv=[1];
        % incoming flux values
        inc(1:sn)   = 0;
    % coarse mesh limit - DSA testing
    case 11
        c = 0.9999;
        % number of elements per zone
        nel_zone = 20;
        % width of each zone
        width_zone = 1;
        % sigt/sigs per zone
        mfp = 100;
        sigt=mfp/(width_zone/nel_zone);
        sigs=c*sigt;
        % volumetric source value, per zone
        qv=[1];
        % incoming flux values
        inc(1:sn)   = 0;
    otherwise
        error('case ID unknown');

end

% sanity checks
if length(nel_zone) ~= length(width_zone)
    error('length(nel_zone) ~= length(width_zone)');
end
if length(nel_zone) ~= length(qv)
    error('length(nel_zone) ~= length(qv)');
end
if length(nel_zone) ~= length(sigt)
    error('length(nel_zone) ~= length(sigt)');
end
if length(nel_zone) ~= length(sigs)
    error('length(nel_zone) ~= length(sigs)');
end
% save data
dat.sigt = sigt;
dat.D = 1./(3.*sigt);
dat.sigs = sigs;
dat.siga = sigt - sigs;
dat.qv =qv;

% a simple check
aux = sigs./sigt;
i=find(aux>1);
if(~isempty(i))
    error('sigs cannot be > than sigt \nA problem occured with zones %i',i);
end

% total number of elements
nel = sum(nel_zone);
% number of spatial dofs
ndof=nel*(porder+1);

% delta x and iel2zon
tmp_width = [ 0 width_zone];
x=[];
iel2zon=[];
for z=1:length(nel_zone)
    x_zone = linspace(tmp_width(z),tmp_width(z+1),nel_zone(z)+1);
    if isempty(x)
        x = x_zone;
        iel2zon=z*ones(nel_zone(z),1);
    else
        x = [x x_zone(2:end)];
        iel2zon =[ iel2zon; z*ones(nel_zone(z),1)];
    end
end
dx=diff(x);
% save data
npar.porder=porder;
npar.ndof=ndof;
npar.nel=nel;
npar.x=x';
npar.dx=dx';
npar.iel2zon=iel2zon;

% angular surfacic source
qsa = zeros(ndof*sn,1);
% add the mu vales to the surfacic source
for k=sn/2+1:sn % positive dir
    qsa(ndof*(k-1)+1)=inc(k);
end
for k=1:sn/2 % negative dir
    qsa(ndof*k)=inc(k);
end
z=kron(abs(mu)',ones(ndof,1));
qsa=qsa.*z;

% element dofs
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=gn(iel-1,:) + npar.porder+1 ;
end
npar.gn = gn;
mats.global.qsa = qsa;