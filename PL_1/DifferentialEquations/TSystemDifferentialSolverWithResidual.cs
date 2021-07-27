// Решение/работа с системой дифференциальных уравнений с n переменными и n неизвестными 
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Решение/работа с системой дифференциальных уравнений с n переменными и n неизвестными 
    /// </summary>
    public partial class TDifferentialSolver
    {
        //------------------------------------------------------------
        /// <summary>
        /// Решение системы дифференциальных уравнений с n переменными и n неизвестными  методом Рунге-Кутты 4ого порядка 
        /// </summary>
        /// <param name="Equation">Решаемая система уравнений и настройки решателя</param>
        /// <returns>Результат решения</returns>
        public static TSystemResultDifferential SolveSystemResidualFourRungeKutta(ISystemDifferentialEquation Equation)
        {
            double X = Equation.Min_X;                              // Крайняя левая точка диапазона "х" 
            double h = Equation.Step;                               // Шаг сетки "h" 
            int t = Equation.Rounding;                              // Округление до нужного знака, после запятой 
            int NumberOfIterations = Equation.CountIterations;      // Количество итераций  ХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХ
            int NumberOfEquations = Equation.CountEquations;        // кол-во урвнений
            List<double> InitArray = Equation.InitArray;            // Начальные значения Y для системы
            TSystemResultDifferential ResultSystemDifferential = new TSystemResultDifferential();


            //List<double> Xs = new List<double>();
            //List<double[]> Ys = new List<double[]>();

            // Рабочие переменные
            double[] Coefs1 = new double[NumberOfEquations];// число 1-ыx коэф. метода по числу уравнений
            double[] Coefs2 = new double[NumberOfEquations];// число 2-ыx коэф. метода по числу уравнений
            double[] Coefs3 = new double[NumberOfEquations];// число 3-иx коэф. метода по числу уравнений
            double[] Coefs4 = new double[NumberOfEquations];// число 4-ыx коэф. метода по числу уравнений

            double[] Y2 = new double[NumberOfEquations]; // число переменных для 2-го коэф. включая независимую
            double[] Y3 = new double[NumberOfEquations];// число переменных для 3-го коэф. включая независимую
            double[] Y4 = new double[NumberOfEquations];// число переменных для 4-го коэф. включая независимую

            //копируем начальные значения игреков в массив, который будет использоваться для вычислений
            double[] Y = new double[NumberOfEquations];//Массив всех Y
            for (int k = 0; k < Equation.InitArray.Count; k++)
            {
                Y[k] = InitArray[k];
            }

            //Закидываем первые значения в резалт
            TPointSystemDifferential PointSystemDifferentialInitial = new TPointSystemDifferential
            {
                Result = new double[NumberOfEquations],
                IndexIteration = 0,
                X = X
            };
            for (int i = 0; i < Y.Length; i++)
            {
                PointSystemDifferentialInitial.Result[i] = Y[i];
            }
            ResultSystemDifferential.SystemPoints.Add(PointSystemDifferentialInitial);

            //Xs.Add(X);
            //Ys.Add(Y);//добавление массивов значений игреков
            for (int i = 0; i < NumberOfIterations; i++)
            {
                TPointSystemDifferential PointSystemDifferential = new TPointSystemDifferential
                {
                    Result = new double[NumberOfEquations],
                    IndexIteration = i + 1,
                    Coeffs = new List<double[]>()
                };

                //Вычисляем новые значения Y и Coeffs для шага h
                var Result =  CalculateValuesOf_Y(X, Y, h, Equation);

                //прибавляем шаг к Х
                X += h;

                //Записываем новые значения в поинт, а далее поинт в резалт
                PointSystemDifferential.X = X;
                for (int j = 0; j < Y.Length; j++)
                {
                    PointSystemDifferential.Result[j] = Result.Y[j];
                }
                PointSystemDifferential.Coeffs.Add(Result.Coefs1);
                PointSystemDifferential.Coeffs.Add(Result.Coefs2);
                PointSystemDifferential.Coeffs.Add(Result.Coefs3);
                PointSystemDifferential.Coeffs.Add(Result.Coefs4);
                ResultSystemDifferential.SystemPoints.Add(PointSystemDifferential);

                //if (i > 0)      // перепроверить
                //    Xs.Add(X);  // перепроверить
                //Ys.Add(Y);      // перепроверить


                //+++++++++++++++++++++++++++++++++++++++//
                //          ИЗМЕНЕНИЕ ШАГА               //
                //var StepCorrectionCoef = CalcStepCorrectionCoef(Xs, Ys, h);
                //if (StepCorrectionCoef == 0d)
                //    return ResultSystemDifferential;
                //h = h * 0.1d / StepCorrectionCoef;
                //+++++++++++++++++++++++++++++++++++++++//



                
            }
            // Вернуть результат
            return ResultSystemDifferential;
        }
        /// <summary>
        /// Метод вычисления новых значений Y с шагом в h
        /// </summary>
        /// <param name="X">Старое значение Х (для Х + h будут вычисляться Y)</param>
        /// <param name="Y">Старые значения Y (на основании которых будут найдены новые значения Y)</param>
        /// <param name="h">Значение "шага"</param>
        /// <param name="Equation">Решаемая система уравнений и настройки решателя</param>
        /// <returns></returns>
        private static (double[] Y, double[] Coefs1, double[] Coefs2, double[] Coefs3, double[] Coefs4) CalculateValuesOf_Y(double X, double[] Y, double h, ISystemDifferentialEquation Equation)
        {
            int NumberOfEquations = Equation.CountEquations;

            double Kx2_3 = X + h / 2;
            double Kx4 = X + h;

            var Coefs1 = Equation.ComputeEquation(X, Y);

            // Находим значения переменных для второго коэф. 
            var Y2 = Do_some_magic(Y, Coefs1, NumberOfEquations, true);
            var Coefs2 = Equation.ComputeEquation(Kx2_3, Y2);

            // Находим значения переменных для третьго коэф.
            var Y3 = Do_some_magic(Y, Coefs2, NumberOfEquations, true);
            var Coefs3 = Equation.ComputeEquation(Kx2_3, Y3);

            // Находим значения переменных для 4 коэф.
            var Y4 = Do_some_magic(Y, Coefs3, NumberOfEquations, false);
            var Coefs4 = Equation.ComputeEquation(Kx4, Y4);

            // Находим новые значения переменных включая независимую    
            for (int k = 0; k < NumberOfEquations; k++)
            {
                Y[k] += (1.0 / 6.0) * (Coefs1[k] + 2 * (Coefs2[k] + Coefs3[k]) + Coefs4[k]);
            }
            return (Y, Coefs1, Coefs2, Coefs3, Coefs4);
        }


        /// <summary>
        /// ✧･ﾟ: *✧･ﾟ:*Makes precious magic*:･ﾟ✧*:･ﾟ✧
        /// </summary>
        /// <param name="Y"></param>
        /// <param name="Coefs">Массив соответствующих коэффициентов</param>
        /// <param name="NumberOfEquations">Число уравнений. Отправляется для выполнения цикла соответствующее кол-во раз</param>
        /// <param name="This_is_2nd_or_3rd_Y">правда если при вызове название данной булиновой переменной верно</param>
        /// <returns></returns>
        private static double[] Do_some_magic (double[] Y, double[] Coefs, int NumberOfEquations, bool This_is_2nd_or_3rd_Y)
        {
            double[] Y_OUT = new double[NumberOfEquations];
            double K234 = new double(); //коэффициент, который только и отличается при расчете У2/У3/У4
            if (This_is_2nd_or_3rd_Y)   //определяется значение этого коэфа
                K234 = 2d;              
            else
                K234 = 1d;
            for (int i = 0; i < NumberOfEquations; i++) 
            {
                Y_OUT[i] = Y[i] + Coefs[i] / K234;       //выполняется прямое предназначение метода
            }
            return Y_OUT;
        }
        //------------------------------------------------------------
        /// <summary>
        /// Вывести отладочную информацию в консоль и в файл если задано имя
        /// </summary>
        /// <param name="Result">Результат решения системы дифф. уравнений</param>
        /// <param name="FileName">Имя файла</param>
        public static void Debug(TSystemResidualResultDifferential Result, string FileName = "")
        {
            // В файл
            if (FileName.Length > 0) File.WriteAllText(FileName, Result.ToString());
            // В консоль
            Console.WriteLine(Result.ToString());
        }
 //------------------------------------------------------------
        /// <summary>
        /// Простой пример системы дифференциальных уравнений  и ее решения 
        /// <returns>Результат решения</returns>
        public static TSystemResultDifferential Example_dN_Residual()
        {

            // Создаем систему уравнений, которая должна решаться и задаем ее параметры 
            ISystemDifferentialEquation Equation = new TEquation_dN()
            {
                Equation = new AEquation_dN((X, Y) =>
                {
                    double[] FunArray = new double[5];//кол-во элементов в массиве должно быть = кол-во уравнений


                    //---------------------------------------------------------------
                    //задаются уравнения
                    FunArray[0] = (X + Y[0] + Y[1] + Y[2] + Y[3] + Y[4]) ;
                    FunArray[1] = (X + 2 * Y[0] + Y[1] + Y[2] + Y[3] + Y[4]);
                    FunArray[2] = (X + Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    FunArray[3] = (5 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    FunArray[4] = (2 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //---------------------------------------------------------------

                    return FunArray;   // интегрируемая система
                }), 
                InitArray = new List<double> { 1, 1, 1, 1, 1, },
                CountIterations = 10,
                Min_X = 0,
                Rounding = 3,
                Step = 1,
                CountEquations = 5
            };
            // Решаем
            return SolveSystemResidualFourRungeKutta(Equation);
        }

        public static double CalcStepCorrectionCoef(List<double> Xs, List<double[]> Ys, double Step)
        {
            int iteration = Ys.Count() - 1;
            
            //указываем граничные значения
            List<TResidual> Residuals_Attributes = new List<TResidual>();
            TResidual R1 = new TResidual
            {
                Name = "H",
                Value = 20000d,
                Accuracy = 1d
            };

            TResidual R2 = new TResidual
            {
                Name = "Cx",
                Value = 2d,
                Accuracy = 0.01d
            };
            //составляем лист граничных значений
            Residuals_Attributes.Add(R1);
            Residuals_Attributes.Add(R2);
            
            //задается метод нахождения пераметров по которым мы смотрим завершать ли вычисление
            var Residuals_calculation = new AEquation_dN((X, Y) =>
            {
                double[] Attribute_Arr = new double[Residuals_Attributes.Count()];//кол-во элементов в массиве должно быть = кол-во методов вычисления
                //---------------------------------------------------------------
                //задаются уравнения по которым будут находиться невязки
                Attribute_Arr[0] = X + Y[0] + Y[4];
                Attribute_Arr[1] = Y[1] + Y[2] + Y[3] + Y[4];
                //---------------------------------------------------------------
                return Attribute_Arr;  //возвращается массив с вычисленными парметрами
            });

            //вычисление текущих значиений по которым идет вычисление невязок для двух различных итераций
            var ValuesOfAttributes_1 = Residuals_calculation(Xs[iteration], Ys[iteration]); //Вычисляются значения параметров для текущей итерации
            var ValuesOfAttributes_0 = Residuals_calculation(Xs[iteration - 1], Ys[iteration - 1]);//Вычисляются значения параметров для преыдущей итерации
            
            //объявление и вычисление изменения невязкок за последнюю итерацию и значения самой невязки
            var Residuals_delta = new double[Residuals_Attributes.Count()];// массив значений изменения невязки
            var Residuals = new double[Residuals_Attributes.Count()];//массив значений истинных невязок 
            for (int i = 0; i < Residuals_delta.Length; i++)
            {
                Residuals_delta[i] = ValuesOfAttributes_0[i] - ValuesOfAttributes_1[i];
                Residuals[i] = Residuals_Attributes[i].Value - ValuesOfAttributes_1[i];
            }

            //проверка, не вошло ли уже значение в область или уже перешагнуло через нее
            for (int i = 0; i < Residuals_Attributes.Count(); i++)
            {
                if ((Math.Abs(Residuals[i]) <= Residuals_Attributes[i].Accuracy) || (ValuesOfAttributes_1[i] > Residuals_Attributes[i].Value + Residuals_Attributes[i].Accuracy))//по хорошему бы переработать для случая, когда условие выхода отрицательное или когда начальное значение параметра больше чем условие выхода и идет спуск к выходному
                    return 0d;
            }

            //объявление массива для корректирующих коэффициентов & заполнение массива корректирующих коэф
            double[] Correction_Coefficients = new double[Residuals_Attributes.Count()];
            for (int i = 0; i < Correction_Coefficients.Count(); i++)
            {
                Correction_Coefficients[i] = Residuals_delta[i] / Residuals_Attributes[i].Accuracy;
            }

            //поиск максимального корректирующего коэф
            double CC_max = new double();
            CC_max = Correction_Coefficients[0];
            for (int i = 0; i < Correction_Coefficients.Length; i++)
            {
                if (Correction_Coefficients[i] > CC_max)
                    CC_max = Correction_Coefficients[i];
            }
            return CC_max;
        }
        //-----------------------------------------------------------
    }
}


//сделать цикл выч У и подбор шага.
//откорректить вывод, когда метод нахождения коэффициента выводит 0
//сделать начало коррекции шага не сразу, а прямо перед завершением.