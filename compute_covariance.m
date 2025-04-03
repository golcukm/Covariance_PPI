function cov_wtb = compute_covariance(data)
    m_data = mean(data);
    data_b = data - m_data;

    x_data = data_b(:,1:3:end);
    y_data = data_b(:,2:3:end);
    z_data = data_b(:,3:3:end);

    cov_x_data = cov(x_data);
    cov_y_data = cov(y_data);
    cov_z_data = cov(z_data);

    cov_data = cov_x_data + cov_y_data + cov_z_data;
    cov_wtb = cov_data;
    for k = 1:size(cov_data, 1)
        cov_wtb(k, k) = 0;
    end
end

