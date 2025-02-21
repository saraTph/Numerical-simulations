function NL = Nonlinearity2BodyInt(g11, g22, g12)
NL = cell(2);
NL{1,1} = @(Phi,X,Y) (g11*(abs(Phi{1})).^2 + g12*(abs(Phi{2})).^2);
NL{1,2} = @(Phi,X,Y) 0;
NL{2,1} = @(Phi,X,Y) 0;
NL{2,2} = @(Phi,X,Y) (g22*(abs(Phi{2})).^2 + g12*(abs(Phi{1})).^2);
end