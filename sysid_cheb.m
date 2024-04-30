function [Phat, Qhat, rad] = sysid_cheb (Yp, Ym, Um, Phi, p, m, L, M, Norm,k)

%======================================================================================================

% Help:
%
% Fuction inputs:   |   Yp              output data matrix Y_+
%                   |   Ym              output data matrix Y_-
%                   |   Um              input  data matrix U_-
%                   |   Phi             Noise model
%                   |   p               output dimension
%                   |   m               input  dimension
%                   |   L               system lag
%                   |   M               M<=L
%                   |   Norm            given norm (choose among 'Frobenius', 'spectral', 'nuclear', 'Ky Fan', 'Schatten'.
%                   |   k               If 'Ky Fan' or 'Schatten' is chosen for Norm, then this input is needed to specify which
%                                       KyFan k-norm or Schatten p-norm is intended. Otherwise, leave it empty.
%
%
% Function outputs: |   Phat            Chebyshev center -- estimation for P_true
%                   |   Qhat            Chebyshev center -- estimation for Q_true
%                   |   rad             Chebyshev radius -- with respect to the specified norm

%======================================================================================================


D=[eye(p)            Yp
   zeros(L*p,p)     -Ym
   zeros((M+1)*m,p) -Um];

N=D*Phi*D';

N11=N(1:p,1:p);
N22=N(p+1:end,p+1:end);
N12=N(1:p,p+1:end);

N_N22=N11-N12*(N22\N12');

cent=-N12/N22;

Phat=cent(:,1:L*p);
Qhat=cent(:,L*p+1:end);

sigma=svd(sqrt(inv(-N22))).*[svd(sqrt(N_N22));0;0];

switch nargin
    case 9

        if strcmp(Norm,'Frobenius') || strcmp(Norm,'frobenius')

            rad=norm(sigma);

        elseif strcmp(Norm,'spectral')

            rad=max(sigma);

        elseif strcmp(Norm,'nuclear')

            rad=sum(sigma);

        else

            msg = "Choose the norm among the following options: 'Frobenius', 'spectral', 'nuclear', 'Ky Fan', 'Schatten'. If 'Ky Fan' or 'Schatten' is chosen, then one more input is needed to specify which KyFan k-norm or Schatten p-norm is intended.";
            error(msg)

        end

    case 10

        if k>p

            msg = "For KyFan k-norms (or Schatten p-norms), k (or p) must be less than or equal to the dimension of the output y.";
            error(msg)

        end

        if strcmp(Norm,'KyFan') || strcmp(Norm,'kyfan') || strcmp(Norm,'Ky Fan') || strcmp(Norm,'ky fan') || strcmp(Norm,'Kyfan') || strcmp(Norm,'kyFan') || strcmp(Norm,'Ky fan') || strcmp(Norm,'ky Fan')

            rad=0;

            for i=1:k

                rad=sigma(i)+rad;

            end

        elseif strcmp(Norm,'Schatten') || strcmp(Norm,'schatten')

            rad2=0;

            for i=1:k

                rad2=sigma(i)^2+rad2;

            end

            rad=sqrt(rad2);

        else

            msg = "Choose the norm among the following options: 'Frobenius', 'spectral', 'nuclear', 'Ky Fan', 'Schatten'. If 'Ky Fan' or 'Schatten' is chosen, then one more input is needed to specify which KyFan k-norm or Schatten p-norm is intended.";
            error(msg)

        end

    otherwise

        msg = "Something is wrong with the specified norm. See the help section.";
        error(msg)

end


