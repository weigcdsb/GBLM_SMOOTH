% Point-process adaptive smoothing w/ Poisson likelihood (log-link)
%  filtering via Eden et al. Neural Comp 2004
%  then a backward pass based on Rauch-Tung-Striebel

function [b,W,lam] = ppafilt_poissexp(n,X,b0,W0,F,Q,offset)

if nargin<7, offset=n*0; end
dt = 1e-3;

% Preallocate
b   = zeros(length(b0),length(n));
W   = zeros([size(W0) length(n)]);
lam = n*0;

% Initialize
b(:,1)   = b0;
W(:,:,1) = W0;
lam(1)   = exp(X(1,:)*b0 + offset(1));

bpred = b;
Wpred = W;

I = eye(size(X,2));

% Forward-Pass (Filtering)
for i=2:length(n)
    bpred(:,i) = F*b(:,i-1);
    lam(i) = exp(X(i,:)*bpred(:,i) + offset(i));
    Wpred(:,:,i) = F*W(:,:,i-1)*F' + Q;
    
    Wpostinv = inv(Wpred(:,:,i)) + X(i,:)'*(lam(i)*dt)*X(i,:);
    W(:,:,i) = inv(Wpostinv);
    
    b(:,i)  = bpred(:,i) + W(:,:,i)*X(i,:)'*(n(i)-lam(i)*dt);
    
    [~, msgid] = lastwarn;
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        return;
    end
end

% Wfilt = W;
% 
% % Backward-Pass (RTS)
% for i=(length(n)-2):-1:1
%     Wi = inv(Wpred(:,:,i+1));
%     Fsquig = inv(F)*(I-Q*Wi);
%     Ksquig = inv(F)*Q*Wi;
%     
%     b(:,i)=Fsquig*b(:,i+1) + Ksquig*bpred(:,i+1);
%     C = W(:,:,i)*F'*Wi;
%     W(:,:,i) = W(:,:,i) + C*(W(:,:,i+1)-Wpred(:,:,i+1))*C';
% end
% 
% A1=zeros(length(b0));
% A2=zeros(length(b0));
% Eznzn_prev = W(:,:,1) + b(:,1)*b(:,1)';
% for i=2:length(n)
%     Eznzn1 = Wfilt(:,:,i-1)*F*inv(F*Wfilt(:,:,i-1)*F'+Q)*W(:,:,i) + b(:,i)*b(:,i-1)';
%     Eznzn = W(:,:,i) + b(:,i)*b(:,i)';
%     A1 = A1+Eznzn1;
%     A2 = A2+Eznzn_prev;
%     Eznzn_prev = Eznzn;
% end
% Anew=A1*inv(A2);
% 
% % Anew=F; % ignore update for linear term
% Qnew=zeros(length(b0));
% Eznzn_prev = W(:,:,1) + b(:,1)*b(:,1)';
% for i=2:length(n)
%     Eznzn1 = Wfilt(:,:,i-1)*F*inv(F*Wfilt(:,:,i-1)*F'+Q)*W(:,:,i) + b(:,i)*b(:,i-1)';
%     Ezn1zn = Wfilt(:,:,i-1)*F*inv(F*Wfilt(:,:,i-1)*F'+Q)*W(:,:,i) + b(:,i-1)*b(:,i)';
%     Eznzn = W(:,:,i) + b(:,i)*b(:,i)';
%     Qnew = Qnew + Eznzn - Anew*Ezn1zn' - Eznzn1*Anew' + Anew*Eznzn_prev*Anew';
%     Eznzn_prev = Eznzn;
% end
% Qnew=Qnew/(length(n)-1);
% keyboard
