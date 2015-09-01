% sym2fun converts a symbolic expression to a function handle.
%
% USAGE:
% ======
% [fun] = sym2fun(expr,var1,var2,...)
%
% INPUTS:
% =======
% expr ... symbolic expression (potentially a vector)
% var1 ... vector containing symbolic variables
% var2 ... vector containing symbolic variables
% ...
%
% Outputs:
% ========
% fun ... function handle which provides the evaluation oe expr given
%           numerical values for var1, ...
%
% 2008/09/23 Steffen Waldherr
% 2012/05/16 Jan Hasenauer

function [fun] = sym2fun(varargin)

%% CHECK AND ASSIGN INPUTS
if nargin >= 2
    expr = varargin{1};
    % Number of variable vectors
    nvar = nargin-1;
else
    error('sym2fun requires at least two input arguments.');
end

%% INITIALIZATION
exp_var = 'theta';

try
%% (1) NO VARIABLE
% This portion of the code handles cases in which the expression expr
% does not contain any variables.

% Check of expr can be converted to numerical expression
%  yes -> no variables contained
%  no  -> variables contained (an error is produced an the algorithm 
%           jumps to the next part (2))

% Construct string expression
expr_num = double(expr);
exprstr = ['@(' exp_var ') ['];
for i = 1:length(expr_num)
    exprstr = [exprstr ';' num2str(expr_num(i))];
end
exprstr = [exprstr ']'];

catch
%% (1) VARIABLES

% Construction of variable vectors and substitution in expression
for k=1:nvar
	x{k} = sym([]);
	argsym = varargin{k+1};
	if isempty(argsym)
        x{k} = [];
    end
	for i=1:size(argsym,1)
        for j=1:size(argsym,2)
        	if nvar == 1
                x{k}(i,j) = sym([exp_var '_' num2str(i) '_' num2str(j)]);
            else
                x{k}(i,j) = sym([exp_var num2str(k) '_' num2str(i) '_' num2str(j)]);
            end
            if ~isreal(argsym(i,j))
                expr = subs(expr,argsym(i,j),x{k}(i,j));
            end
        end
    end
end

% Construction of string from expression
exprstr = '@(?) [';
for d1=1:size(expr,1)
	for d2=1:size(expr,2)
		next = char(expr(d1,d2));
		for k=1:nargin-1
        	for i=size(x{k},1):-1:1
                for j=size(x{k},2):-1:1
                    if nvar == 1
                        if size(x{k},2) == 1
                            next = strrep(next,char(x{k}(i,j)),['?' '(' num2str(i) ')']);
                        else
                            next = strrep(next,char(x{k}(i,j)),['?' '(' num2str(i) ',' num2str(j) ')']);
                        end
                    else
                        if size(x{k},2) == 1
                            next = strrep(next,char(x{k}(i,j)),['?' num2str(k) '(' num2str(i) ')']);
                        else
                            next = strrep(next,char(x{k}(i,j)),['?' num2str(k) '(' num2str(i) ',' num2str(j) ')']);
                        end
                    end
                end
            end
        end
		exprstr = [exprstr next];
		if d2<size(expr,2)
			exprstr = [exprstr ','];
		elseif d1<size(expr,1)
			exprstr = [exprstr ';'];
        end
    end
end
exprstr = [exprstr '];'];
for k=1:nvar
    if nvar == 1
        exprstr = strrep(exprstr,['?'],[exp_var]);
    else
        exprstr = strrep(exprstr,['?' num2str(k)],[exp_var num2str(k)]);
    end
end
if nvar == 1
    argstr = [exp_var];
else
    argstr = [exp_var '1'];
end
for k=2:nvar
	argstr = [argstr ', ' exp_var num2str(k)];
end
exprstr = strrep(exprstr,'?0',argstr);

end

% Compile string expression to function
fun = eval(exprstr);

