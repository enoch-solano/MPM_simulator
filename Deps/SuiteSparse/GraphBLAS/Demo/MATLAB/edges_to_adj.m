function C = edges_to_adj (E)
%EDGES_TO_ADJ convert incidence matrix to adjacency matrix
%
% usage: C = edges_to_adj (E)
%
% Each column of E is an edge (i,j) where E([i j],:) = 1,
% Each column must have exactly two or zero entries.
%
% C is a symmetric binary matrix with no self edges.

% check input E
[i j x] = find (E) ;
if (any (x ~= 1))
    error ('invalid incidence matrix: not binary') ;
end
d = full (sum (E)) ;
if (~all ((d == 2) | (d == 0)))
    error ('invalid incidence matrix: column degree') ;
end

% convert each edge into a clique of size 2, and sum up the cliques
C = spones (E*E') ;

% remove the diagonal
C = C - diag (diag (C)) ;
% C = tril (C,-1) + triu (C,1) ;

