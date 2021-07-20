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
        public static TSystemResidualResultDifferential SolveSystemResidualFourRungeKutta(ISystemDifferentialEquation Equation)
        {
            double X = Equation.Min_X;                              // Крайняя левая точка диапазона "х" 
            double h = Equation.Step;                               // Шаг сетки "h" 
            int t = Equation.Rounding;                              // Округление до нужного знака, после запятой 
            int CountOfIterations = Equation.CountIterations;       // Количество итераций  ХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХХ
            int CountEquation = Equation.CountEquatiuon;            // кол-во урвнений
            List<double> InitArray = Equation.InitArray;            // Начальные значения Y для системы
            TSystemResidualResultDifferential ResultSystemDifferential = new TSystemResidualResultDifferential();


            List<double> Xs = new List<double>();
            List<double[]> Ys = new List<double[]>();

            // Рабочие переменные
            double[] Coefs1 = new double[CountEquation];// число 1-ыx коэф. метода по числу уравнений
            double[] Coefs2 = new double[CountEquation];// число 2-ыx коэф. метода по числу уравнений
            double[] Coefs3 = new double[CountEquation];// число 3-иx коэф. метода по числу уравнений
            double[] Coefs4 = new double[CountEquation];// число 4-ыx коэф. метода по числу уравнений

            double[] Y2 = new double[CountEquation]; // число переменных для 2-го коэф. включая независимую
            double[] Y3 = new double[CountEquation];// число переменных для 3-го коэф. включая независимую
            double[] Y4 = new double[CountEquation];// число переменных для 4-го коэф. включая независимую

            //копируем начальные значения игреков в массив, который будет использоваться для вычислений
            double[] Y = new double [CountEquation];//Массив всех Y
            for (int k = 0; k < Equation.InitArray.Count; k++)
            {
                Y[k] = InitArray[k];
            }

            Xs.Add(X);
            Ys.Add(Y);//добавление массивов значений игреков
            for (int i = 0; i < CountOfIterations; i++)
            {
                TPointSystemResidualDifferential PointSystemDifferential = new TPointSystemResidualDifferential
                {
                    Result = new double[CountEquation],
                    IndexIteration = i,
                    Coeffs = new List<double[]>()
                };

                Coefs1 = Equation.ComputeEquation(X, Y);
                // Находим значения переменных для второго коэф.    
                double Kx2_3 = X + h / 2;
                double Kx4 = X + h;


                for (int k = 0; k < CountEquation; k++)
                {
                    Y2[k] = Y[k] + Coefs1[k] / 2;
                }
                Coefs2 = Equation.ComputeEquation(Kx2_3, Y2);

                // Находим значения переменных для третьго коэф.
                for (int k = 0; k < CountEquation; k++)
                {
                    Y3[k] = Y[k] + Coefs2[k] / 2;
                }

                Coefs3 = Equation.ComputeEquation(Kx2_3, Y3);

                // Находим значения переменных для 4 коэф.

                for (int k = 0; k < CountEquation; k++)
                {
                    Y4[k] = Y[k] + Coefs3[k];
                }
                Coefs4 = Equation.ComputeEquation(Kx4, Y4);

                // Находим новые значения переменных включая независимую    

                for (int k = 0; k < CountEquation; k++)
                {
                    Y[k] += (1.0 / 6.0) * (Coefs1[k] + 2 * (Coefs2[k] + Coefs3[k]) + Coefs4[k]);
                }
                PointSystemDifferential.X = X;
                // Результат  иттерации:
                for (int j = 0; j < CountEquation; j++)
                {
                    PointSystemDifferential.Result[j] = Y[j];
                }
                PointSystemDifferential.Coeffs.Add(Coefs1);
                PointSystemDifferential.Coeffs.Add(Coefs2);
                PointSystemDifferential.Coeffs.Add(Coefs3);
                PointSystemDifferential.Coeffs.Add(Coefs4);
                ResultSystemDifferential.SystemPoints.Add(PointSystemDifferential);

                if (i > 0)
                    Xs.Add(X);
                Ys.Add(Y);
                h = ReCalcStep(Xs, Ys, h);
                if (h == 0d)
                    return ResultSystemDifferential;



                X += h;
            }
            // Вернуть результат
            return ResultSystemDifferential;
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
        public static TSystemResidualResultDifferential Example_dN_Residual()
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

                    return FunArray;   // интегрируемая система + вычесленные невязки
                }), 
                InitArray = new List<double> { 1, 1, 1, 1, 1, },
                CountIterations = 10,
                Min_X = 0,
                Rounding = 3,
                Step = 1,
                CountEquatiuon = 5
            };
            // Решаем
            return SolveSystemResidualFourRungeKutta(Equation);
        }

        public static double ReCalcStep(List<double> Xs, List<double[]> Ys, double Step)
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
            
            //zadayetsya metod vichislenie nevyazok
            var Residuals_calculation = new AEquation_dN((X, Y) =>
            {
                double[] Residuals_arr = new double[Residuals_Attributes.Count()];//кол-во элементов в массиве должно быть = кол-во невязок
                //---------------------------------------------------------------
                //задаются уравнения
                Residuals_arr[0] = X + Y[0] + Y[4];
                Residuals_arr[1] = Y[1] + Y[2] + Y[3] + Y[4];
                //---------------------------------------------------------------
                return Residuals_arr;   // интегрируемая система + вычесленные невязки
            });

            var Residuals_1 = Residuals_calculation(Xs[iteration], Ys[iteration]);
            var Residuals_0 = Residuals_calculation(Xs[iteration - 1], Ys[iteration - 1]);
            var Residual_delta = new double[Residuals_Attributes.Count()];
            var Residual_Rmng = new double[Residuals_Attributes.Count()];
            for (int i = 0; i < Residual_delta.Length; i++)
            {
                Residual_delta[i] = Residuals_0[i] - Residuals_1[i];
                Residual_Rmng[i] = Residuals_Attributes[i].Value - Residuals_1[i];
            }

            //объявление массива для корректирующих коэффициентов & заполнение массива корректирующих коэф
            double[] Correction_Coefficients = new double[Residuals_Attributes.Count()];
            for (int i = 0; i < Correction_Coefficients.Count(); i++)
            {
                Correction_Coefficients[i] = Residual_delta[i] / Residuals_Attributes[i].Accuracy;
            }

            double CC_max = new double();
            CC_max = Correction_Coefficients[0];
            for (int i = 0; i < Correction_Coefficients.Length; i++)
            {
                if (Correction_Coefficients[i] > CC_max)
                    CC_max = Correction_Coefficients[i];
            }


            //проверка, не вошло ли уже значение в область
            for (int i = 0; i < Residuals_Attributes.Count(); i++)
            {
                if ((Math.Abs(Residuals_1[i]) < Residuals_Attributes[i].Accuracy) || (Residuals_1[i] + Residuals_Attributes[i].Accuracy > Residuals_Attributes[i].Value))//pereproverit'
                    return 0d;
            }
            return CC_max;
        }
        //-----------------------------------------------------------
    }
}



//в солвер надо добавить выход если шаг = 0
//надо реализовать изменение шага