function A = get_A(theta)
    theta_sub = theta(1:end-1);
    theta_1 = theta_sub(1:length(theta_sub)/2);
    theta_2 = theta_sub(length(theta_sub)/2+1:end);
    A = [theta_1(1) - theta_1(2), theta_1(3);theta_2(3),theta_2(1)-theta_2(2)];
end