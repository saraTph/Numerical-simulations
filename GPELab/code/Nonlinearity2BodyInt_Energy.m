function NLE = Nonlinearity2BodyInt_Energy(g11, g22, g12)
NLE = cell(2);
NLE{1,1} = @(Phi,X,Y) (1/2)*(g11*(abs(Phi{1})).^2 + g12*(abs(Phi{2})).^2);
NLE{1,2} = @(Phi,X,Y) 0;
NLE{2,1} = @(Phi,X,Y) 0;
NLE{2,2} = @(Phi,X,Y) (1/2)*(g22*(abs(Phi{2})).^2 + g12*(abs(Phi{1})).^2);
end