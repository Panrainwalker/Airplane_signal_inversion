% -------------------------------------------------------------------------
% update for right hand side (rhs) derivatives; backward
%
%
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: May 31, 2025
% -------------------------------------------------------------------------
%%%%%%

% hW(src_iz,src_ix,1) = hW(src_iz,src_ix,1) + stf(it) * 1e6;
% hW(src_iz,src_ix,2) = hW(src_iz,src_ix,2) + stf(it) * 1e6;

% hW(src_iz,src_ix,3) = hW(src_iz,src_ix,3) - stf(it)*1e6;
% hW(M.izs,M.ixs,4) = hW(M.izs,M.ixs,4) - stf(it)*M.mzz;
% hW(src_iz,src_ix,5) = hW(src_iz,src_ix,5) - stf(it)*1e6;


% Loop through all sources
ix = round(src_ix_start + v0 * it * dt / dx);

% 边界检查，避免越界
if ix > ndamp && ix < nx - ndamp
    W(src_iz, ix, 3) = W(src_iz, ix, 3) + stf(it) * 1e10;
end