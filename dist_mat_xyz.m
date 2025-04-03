function [dist_mat, dd] = dist_mat_xyz(m_kif1a_k1,coff)
    mx = m_kif1a_k1(1,1:3:end);
    my = m_kif1a_k1(1,2:3:end);
    mz = m_kif1a_k1(1,3:3:end);

    dist_mat = zeros(size(mx,2),size(mx,2));
    dd = zeros(size(mx,2),size(mx,2));
    for i = 1:size(dist_mat,2)
        for j = 1:size(dist_mat,2)
            tmp_dist = sqrt((mx(j) - mx(i))^2 + (my(j) - my(i))^2 + (mz(j) - mz(i))^2);
            dd(i,j) = tmp_dist;
            if tmp_dist > coff
                dist_mat(i,j) = 0;
            else
                dist_mat(i,j) = 1;
            end
        end
    end
end