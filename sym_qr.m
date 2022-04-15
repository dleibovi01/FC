function [Q,R] = sym_qr(A)

[M,N] = size(A);

if M ~= N
    error 'sym_qr only works on square matrices!'
end

Q = sym( zeros( size(A) ) );
R = sym( zeros( size(A) ) );

Q(:,1) = A(:,1);
m = sqrt(Q(:,1).'*Q(:,1));
R(1,1) = m;
Q(:,1) = Q(:,1)/m;

for k=2:N
   Q(:,k) = A(:,k);
   for l = 1:k-1
       R(l,k) = Q(:,l).'*Q(:,k);
   end
   Q(:,k) = Q(:,k) - Q*R(:,k);
   m = sqrt(Q(:,k).'*Q(:,k));
   Q(:,k) = Q(:,k)/m;
   R(k,k) = m;
end