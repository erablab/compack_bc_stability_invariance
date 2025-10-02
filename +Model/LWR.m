classdef LWR < Model.ModelBase	
	properties (SetAccess = private)
		name = 'LWR'
		nelem = 1
	end
	
	
	methods
		function ret = f(obj, U, d)
			ret = U.*(1 - U);
		end
	end
	
	
	methods(Access=protected)
		function ret = maxEigRect(obj, U, d)
			ret = abs(U);
		end
	end
end