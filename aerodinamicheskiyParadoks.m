classdef AerodynamicParadox < handle                                
                                                                           
% Данный класс относится к моделированию т.н. "аэродинамического парадокса 
% спутника" — попадая в верхние слои атмосферы, космический аппарат, 
% испытывая торможение в разреженном газе, увеличивает при этом скорость 
% своего движения.
% https://kvant.mccme.ru/pdf/2011/04/Travin.pdf
% 
    
    properties 
        
        % Мировые константы                                                
        GM = 6.67.*10.^(-11).*6.*10.^24;   % Гравитационная постоянная
        Rearth = 6400000;   % Радиус Земли

        % Параметры спутника                                               
        SfM = 1; % отн.ед., Площадь миделя
        Hic = 500000; % м, Исходное расстояние от поверхности земли до спутника 
        Fic = 0; % рад, Полярный угол положения спутника в начале отсчета
        dFdtic  % рад/сек, Изменение полярного угла

        % Параметры окружающей среды                                       
        Cx = 2;    % Коэффициент лобового сопротивления
        rho0 = 5.*10.^(-8); % кг/м3, Плотность атмосферы для соответствующей высоты Hic
        H = 40000;   % Высота однородной атмосферы

        % Константы задачи
        tau = 0.1;   % Шаг по сетке времени

    end
    
    methods
        function obj = AerodynamicParadox(name, value)              
                                                                           
            % Конструктор
            arguments (Repeating)                                          
                name {mustBeNoempty, mustBeText}
                value {mustBeNoempty, mustBeNumeric}
            end

            obj.dFdtic = sqrt(obj.GM. / (obj.Hic + obj.Rearth).^3);          
            
            for idx = 1:numel(value)                                         
                obj.(name{idx}) = value{idx};
            end
        clearvars
        end
                            
        function [figureHandle] = plot_velocity_and_trajectory(obj)        
            % Создание графиков для вычисления скорости и траектории движения спутника по орбите

            [radius, angle, radiisTimeDerivative, anglesTimeDerivative] = obj.solve_diff_equation(10); 
            
            figureHandle = figure;
            tileHandle = tiledlayout(1, 2);

            ax1 = nexttile;                                                
            ax2 = nexttile;                                                
            
            % Создание графика для вычисления скорости через изменения угла и радиуса по времени                                                              
            h1 = plot(ax1, sqrt(anglesTimeDerivative.^2 .*radius.^2 ...     
                                + radiisTimeDerivative.^2));      
            
            % Создание графика для наглядности траектории движения спутника в сравнении с окружностью орбиты 
            axes(ax2)
            h2 = polarplot( [0:0.01:2.*pi], ...                             
                           obj.Rearth .*ones(size([0:0.01:2.*pi])), 'r');   
            ax2 = gca;
            hold(ax2, 'on')
            h3 = polarplot(ax2, angle, radius, 'b')
            end
            

        function [R, F, dRdt, dFdt] = solve_diff_equation(obj, periodsAmount)
            % Решение разностной схемы для движения спутника по орбите с трением. 

            timesteps = [0:obj.tau:periodsAmount.*2.*pi./obj.dFdtic]';     

            R = zeros(numel(timesteps), 1);    % Расстояние от центра Земли до спутника
            F = zeros(numel(timesteps), 1);    % Полярный угол положения спутника в инерциальной СК с центром в центре Земли
            dRdt = zeros(numel(timesteps), 1);   % Производная радиуса по времени
            dFdt = zeros(numel(timesteps), 1);   % Производная полярного угла по времени
            
            % 0-ой узел
            R(1) = obj.Hic + obj.Rearth;
            F(1) = obj.Fic;

            % 1-ый узел
            R(2) = obj.Hic + obj.Rearth;
            F(2) = obj.Fic + obj.dFdtic.*obj.tau;

            for k = 2:numel(timesteps)-1                                   
                % 1-ые производные и угол
                dValdt_func = @(Val2, Val1) (Val2 - Val1)./obj.tau;              
                beta_func = @(dRdt, dFdt, R) atan(dRdt./(dFdt.*R));        
                
                % Сила трения 
                prefactor = obj.Cx.*obj.SfM.*obj.rho0./2;               
                Ftr_func = @(dRdt, dFdt, R) prefactor.*(dRdt.^2 + dFdt.^2.*R.^2) .* exp(- (R - obj.Rearth)./obj.H); 
                
                % Производные
                dRdt(k) = dValdt_func(R(k), R(k-1));
                dFdt(k) = dValdt_func(F(k), F(k-1));
                    
                Ftr = Ftr_func(dRdt(k), dFdt(k), R(k));
            
                angleBeta = beta_func(dRdt(k), dFdt(k), R(k));
                
                % Координаты
                R(k+1) = 2.*R(k) - R(k-1) + obj.tau.^2 .*(R(k).*dFdt(k).^2 - obj.GM./R(k).^2 - Ftr.*sin(angleBeta));
                F(k+1) = 2.*F(k) - F(k-1) + obj.tau.^2 ./R(k) .* (-Ftr.*cos(angleBeta) - 2.*dFdt(k) .*dRdt(k));
            end
            
            % Производные последнего узла
            dRdt(k+1) = dValdt_func(R(k+1), R(k));
            dFdt(k+1) = dValdt_func(F(k+1), F(k));
            
            % Удалить 0-ой узел
            R(1) = [];
            F(1) = [];
            dRdt(1) = [];
            dFdt(1) = [];

        end
    end

end












































