classdef ControlledStabilityAndInvarianceSeparate
    %PERIODIC   Periodic boundary condition


    methods (Static)
        function U = updateBoundary(U, mesh, t)
            if mod(t, 1) <= 1e-1
                disp(t)
            end

            ngc = mesh.ngc;
            nx = mesh.nx;
            ndims = mesh.ndims;
            m = 1:ngc;		% Boundary index

            DX = 1 / 100;
            GAMMA_V = 1e-5;
            GAMMA_H = 1e-1;
            u_star = 1/3 * ones(size(U));

            f = @(u) u .* (1 - u);
            F = @(u) u ^ 2 / 2 - u ^ 3 / 3;
            V = @(u) 1 / 2 * sum((u - u_star) .^ 2) * DX;
            H = @(u) (1/4) ^ 2 - sum(u .^ 2) * DX;

            dV_dt = @(u_a, u_b) (u_a - u_star(1)) * f(u_a) ...
                - (u_b - u_star(end)) * f(u_b) ...
                - F(u_a) ...
                + F(u_b);

            dH_dt = @(u_a, u_b) - 2 * u_a * f(u_a) ...
                + 2 * u_b * f(u_b) ...
                + F(u_a) ...
                - F(u_b);

            dv_constr = @(u_a, z1) dV_dt(u_a, z1) + GAMMA_V * V(U);
            dh_constr = @(u_a, z2) -dH_dt(u_a, z2) - GAMMA_H * H(U);

            function [cineq, ceq] = constr_V(ua_z1_z2)
                cineq = [dv_constr(ua_z1_z2(1), ua_z1_z2(2));
                    dh_constr(ua_z1_z2(1), ua_z1_z2(3))];
                ceq = [];
            end

            function [cineq, ceq] = constr_H(ua_z1_z2)
                cineq = [dv_constr(ua_z1_z2(2), ua_z1_z2(1));
                    dh_constr(ua_z1_z2(3), ua_z1_z2(1))];
                ceq = [];
            end

            cost_V = @(ua_z1_z2) ua_z1_z2(1) ^ 2;

            cost_H = @(ub_z1_z2) ub_z1_z2(1) ^ 2;

            ua_z1_z2_star = fmincon(cost_V, [0; 0; 0], ...
                [], [], [], [], [0; 1 / 4; 1 / 4], [3 / 4; 1; 1], ...
                @constr_V, ...
                optimoptions('fmincon','Display','off','MaxIterations',5));

            ub_z1_z2_star = fmincon(cost_H, [0; 0; 0], ...
                [], [], [], [], [1 / 4; 0; 0], [1; 3 / 4; 3 / 4], ...
                @constr_H, ...
                optimoptions('fmincon','Display','off','MaxIterations',5));

            global data
            data.Vvec = [data.Vvec, V(U)];
            data.Hvec = [data.Hvec, H(U)];
            data.u_a = [data.u_a, ua_z1_z2_star(1)];
            data.u_b = [data.u_b, ub_z1_z2_star(1)];

            if ndims >= 1
                U(:, m, :) = ua_z1_z2_star(1);
                U(:, m+ngc, :) = ua_z1_z2_star(1);
                U(:, nx(1)+m, :) = ub_z1_z2_star(2);
                U(:, nx(1)+m+ngc, :) = ub_z1_z2_star(2);
            end
        end
    end
end

