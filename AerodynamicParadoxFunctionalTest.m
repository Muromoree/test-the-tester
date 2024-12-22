classdef (TestTads = {'functionaltest'})AerodynamicParadoxTest < matlab.unitest.TestCase
methods (TestClassTeardown)
        function teardown (~)
        clearvars;
    end
end

methods
    function RadiusAndAngleValues (TestCase)
        %Тест на получение радиусов и углов
        % Создание объекта класса
            obj = AerodynamicParadoxTest();
            
            periodsAmount = 10;  % Задаем количество периодов
            [R, F, dRdt, dFdt] = obj.solvediffequation(periodsAmount);
            
            % Проверка размера выходных данных
            expectedLinearSize = ceil((periodsAmount  2  pi / obj.dFdtic) / obj.tau) + 1;
            testCase.verifySize(R, [expectedLinearSize, 1], 'Size of R array is incorrect.');
            testCase.verifySize(F, [expectedLinearSize, 1], 'Size of F array is incorrect.');
            testCase.verifySize(dRdt, [expectedLinearSize, 1], 'Size of dRdt array is incorrect.');
            testCase.verifySize(dFdt, [expectedLinearSize, 1], 'Size of dFdt array is incorrect.');
            
            % Проверка начальных значений
            expectedInitialRadius = obj.Hic + obj.Rearth;  
            expectedInitialAngle = 0;  
            
            testCase.verifyEqual(R(1), expectedInitialRadius, 'Initial radius value is incorrect.');
            testCase.verifyEqual(F(1), expectedInitialAngle, 'Initial angle value is incorrect.');
        end
        
    function test_DerivativesComputation(testCase)
        % Тест на проверку вычисления производных
            
        [R, F, dRdt, dFdt] = obj.solvediffequation(periodsAmount);
            
        % Проверяем, что производные определяются корректно (для простоты просто проверяем размеры)
        testCase.verifySize(dRdt, [size(R)], 'Size of dRdt is incorrect.');
        testCase.verifySize(dFdt, [size(F)], 'Size of dFdt is incorrect.');
            
        for k = 2:numel(R)-1
            expected_dRdt = (R(k) - R(k-1)) / obj.tau;  
            expected_dFdt = (F(k) - F(k-1)) / obj.tau;  
                
            testCase.verifyEqual(dRdt(k), expected_dRdt, 'dRdt calculation is incorrect at k=%d.', k);
            testCase.verifyEqual(dFdt(k), expected_dFdt, 'dFdt calculation is incorrect at k=%d.', k);
        end

    function test_CreateFigureWithSubplots(testCase)
        % Тест на проверку отображения графиков    

        % Вызываем метод plotvelocityand_trajectory
        figureHandle = obj.plotvelocityand_trajectory();
            
        % Проверяем, что графическое окно создано
        testCase.verifyTrue(isgraphics(figureHandle), 'Figure was not created.');

        % Проверяем, что в фигуре два подграфика
        ax = findobj(figureHandle, 'Type', 'Axes');  % Находим все подграфики
        testCase.verifyEqual(numel(ax), 2, 'There should be two subplots.');  % Проверяем количество подграфиков
            
        % Проверяем, что оба подграфика непустые
        for i = 1:2
            testCase.verifyNotEmpty(ax(i).Children, 'Subplot %d is empty.', i);
         end
    end
end       
