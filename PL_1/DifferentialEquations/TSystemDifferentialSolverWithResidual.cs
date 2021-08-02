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

            //объявляются списки в которых будут храниться значения, которые были вычесленны. 
            List<double> Xs = new List<double>();
            List<double[]> Ys = new List<double[]>();

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

            //Закидываем первые значения в резалт поинт и далее в резалт
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

            List<double[]> Values = new List<double[]>();
            var TrackedValues = CalculateValuesOfTrackedVariables(X, Y);
            double[] Values_arr = new double[TrackedValues.Count()];
            for (int i = 0; i < TrackedValues.Count(); i++)
            {
                Values_arr[i] = TrackedValues[i].CurrentValue;
            }
            Values.Add(Values_arr);
            
            
            for (int i = 0; i < NumberOfIterations; i++)
            {
                Xs.Add(X);// запись значений
                double[] Ys_arr = new double[Y.Length];
                for (int j = 0; j < Y.Length; j++)
                {
                    Ys_arr[j] = Y[j];
                }
                Ys.Add(Ys_arr);//добавление массивов значений игреков

                (double[] Y, double[] Coefs1, double[] Coefs2, double[] Coefs3, double[] Coefs4) Result;
                bool Incorrect = new bool();
                bool TimeToStop = new bool();
                bool StepChanged = new bool();
                TimeToStop = false;
                StepChanged = false;
                do
                {
                    StepChanged = false;
                    Incorrect = true;
                    //Вычисляем новые значения Y и Coeffs для шага h
                    Result = CalculateValuesOf_Y(X, Y, h, Equation);

                    //просчет контрольных значений и невязок, так же хранит точность, название 
                    TrackedValues = CalculateValuesOfTrackedVariables(X + h, Y);

                    if (i < 1)
                        break;
                    
                    double[] Deltas = new double[TrackedValues.Count()];
                    for (int j = 0; j < TrackedValues.Count(); j++)
                    {
                        Deltas[j] =/* Y[j] -*/TrackedValues[j].CurrentValue - Values[Values.Count() - 1][j];
                    }

                    //проверка на необходимость уменьшения шага
                    for (int j = 0; j < TrackedValues.Count(); j++)
                    {
                        if (Deltas[j] > 0)
                        {
                            if (TrackedValues[j].CurrentValue > TrackedValues[j].ControlValue + TrackedValues[j].Accuracy)
                            {
                                h = h / 2d;
                                StepChanged = true;
                                break;
                            }
                            else
                            {
                                //ok ok
                            }
                        }
                        else if (Deltas[j] < 0)
                        {
                            if (TrackedValues[j].CurrentValue < TrackedValues[j].ControlValue - TrackedValues[j].Accuracy)
                            {
                                h = h / 2d;
                                StepChanged = true;
                                break;
                            }
                            else
                            {
                                //ok ok
                            }
                        }
                    }

                    if (StepChanged)
                    {
                        for (int j = 0; j < Y.Length; j++)
                        {
                            Y[j] = Ys[Ys.Count() - 1][j];
                        }
                    }
                    else
                    {
                        Incorrect = false;
                        //проверка, не вошло ли уже значение в область
                        for (int j = 0; j < TrackedValues.Count(); j++)
                        {
                            if (Math.Abs(TrackedValues[j].Residual) <= TrackedValues[j].Accuracy)//по хорошему бы переработать для случая, когда условие выхода отрицательное или когда начальное значение параметра больше чем условие выхода и идет спуск к выходному
                            {
                                TimeToStop = true;
                                break;
                            }
                        }
                    }
                }
                while (Incorrect);


                Values_arr = new double[TrackedValues.Count()];
                for (int j = 0; j < TrackedValues.Count(); j++)
                {
                    Values_arr[j] = TrackedValues[j].CurrentValue;
                }
                Values.Add(Values_arr);

                //прибавляем шаг к Х
                X += h;

                TPointSystemDifferential PointSystemDifferential = new TPointSystemDifferential // перенести в конец цикла---------------------------------
                {
                    Result = new double[NumberOfEquations],
                    IndexIteration = i + 1,
                    Coeffs = new List<double[]> { Result.Coefs1, Result.Coefs2, Result.Coefs1, Result.Coefs3, Result.Coefs4 }
                };
                //Записываем новые значения в поинт, а далее поинт в резалт
                PointSystemDifferential.X = X;
                for (int j = 0; j < Y.Length; j++)
                {
                    PointSystemDifferential.Result[j] = Result.Y[j];
                }
                ResultSystemDifferential.SystemPoints.Add(PointSystemDifferential);

                if (TimeToStop)
                    return ResultSystemDifferential;





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
            var Y2 = CalcYVolumeForCoefs(Y, Coefs1, NumberOfEquations, h, true);
            var Coefs2 = Equation.ComputeEquation(Kx2_3, Y2);

            // Находим значения переменных для третьго коэф.
            var Y3 = CalcYVolumeForCoefs(Y, Coefs2, NumberOfEquations, h, true);
            var Coefs3 = Equation.ComputeEquation(Kx2_3, Y3);

            // Находим значения переменных для 4 коэф.
            var Y4 = CalcYVolumeForCoefs(Y, Coefs3, NumberOfEquations, h, false);
            var Coefs4 = Equation.ComputeEquation(Kx4, Y4);

            // Находим новые значения переменных включая независимую    
            for (int k = 0; k < NumberOfEquations; k++)
            {
                Y[k] += (Coefs1[k] + 2 * (Coefs2[k] + Coefs3[k]) + Coefs4[k]) * h / 6d;
            }
            return (Y, Coefs1, Coefs2, Coefs3, Coefs4);
        }


        /// <summary>
        /// Вычисляет значение Y для вычисление коэффициентов
        /// </summary>
        /// <param name="Y">Старые значения игреков на основании, которых вычисляются коэф., а затем на основании коэф. вычисляются новые значения игреков</param>
        /// <param name="Coefs">Массив соответствующих коэффициентов</param>
        /// <param name="NumberOfEquations">Число уравнений. Отправляется для выполнения цикла соответствующее кол-во раз</param>
        /// <param name="This_is_2nd_or_3rd_Y">правда если при вызове название данной булиновой переменной верно</param>
        /// <returns></returns>
        private static double[] CalcYVolumeForCoefs(double[] Y, double[] Coefs, int NumberOfEquations, double h, bool This_is_2nd_or_3rd_Y)
        {
            double[] Y_OUT = new double[NumberOfEquations];
            double K234 = new double(); //коэффициент, который только и отличается при расчете У2/У3/У4
            if (This_is_2nd_or_3rd_Y)   //определяется значение этого коэфа
                K234 = 2d;
            else
                K234 = 1d;
            for (int i = 0; i < NumberOfEquations; i++)
            {
                Y_OUT[i] = Y[i] + h * Coefs[i] / K234;       //выполняется прямое предназначение метода
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
                    double[] FunArray = new double[1];//кол-во элементов в массиве должно быть = кол-во уравнений


                    //---------------------------------------------------------------
                    //задаются уравнения
                    FunArray[0] = (X * X - 2 * Y[0]);
                    //FunArray[1] = (X + 2 * Y[0] + Y[1] + Y[2] + Y[3] + Y[4]);
                    //FunArray[2] = (X + Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //FunArray[3] = (5 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //FunArray[4] = (2 * X + 2 * Y[0] + 3 * Y[1] + Y[2] + Y[3] + Y[4]);
                    //---------------------------------------------------------------

                    return FunArray;   // интегрируемая система
                }),
                InitArray = new List<double> { 1/*, 1, 1, 1, 1, */},
                CountIterations = 1000,
                Min_X = 0,
                Rounding = 3,
                Step = 0.1,
                CountEquations = 1
            };
            // Решаем
            return SolveSystemResidualFourRungeKutta(Equation);
        }

        private static List<TResidual> CalculateValuesOfTrackedVariables(double X, double[] Y)
        {
            List<TResidual> Residuals = new List<TResidual>();
            TResidual R1 = new TResidual()
            {
                Name = "",
                ControlValue = 1.99d,
                CurrentValue = Y[0],
                Accuracy = 0.02d
            };

            TResidual R2 = new TResidual()
            {
                Name = "",
                ControlValue = 0.3d,
                CurrentValue = X,
                Accuracy = 0.01d
            };

            //Residuals.Add(R1);
            Residuals.Add(R2);

            return Residuals;
        }

      
        public static double CalcStepCorrectionCoef(List<double> Xs, List<double[]> Ys, double Step)
        {
            int iteration = Ys.Count() - 1;
            
            //указываем граничные значения
            List<TResidual> Residuals_Attributes = new List<TResidual>();
            TResidual R1 = new TResidual
            {
                Name = "H",
                ControlValue = 20000d,
                Accuracy = 1d
            };

            TResidual R2 = new TResidual
            {
                Name = "Cx",
                ControlValue = 2d,
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
                Residuals[i] = Residuals_Attributes[i].ControlValue - ValuesOfAttributes_1[i];
            }

            //проверка, не вошло ли уже значение в область или уже перешагнуло через нее
            for (int i = 0; i < Residuals_Attributes.Count(); i++)
            {
                if ((Math.Abs(Residuals[i]) <= Residuals_Attributes[i].Accuracy) || (ValuesOfAttributes_1[i] > Residuals_Attributes[i].ControlValue + Residuals_Attributes[i].Accuracy))//по хорошему бы переработать для случая, когда условие выхода отрицательное или когда начальное значение параметра больше чем условие выхода и идет спуск к выходному
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


//обойти изменение У в do while
//используется изменение Y без возврата самого Y, требуется перепроверить, можно ли вообше не выводить Y?   CalculateValuesOf_Y