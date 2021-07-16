//Вычисление невязок для дифференциального уравнения
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace StandartHelperLibrary.MathHelper
{
    class TResidual
    {
        public double Residual(ISystemDifferentialEquation Equation)
        {
            List<double> Y = Equation.InitArray;                        // Начальные значения Y 
            List<double> Ypol = new List<double>();                     // значения Y, которые задает пользователь
            double h = Equation.Step;                                   //шаг интегрирования
            double CountOfIterations = Equation.CountIterations;        //Количество итераций


            //проверяем шаг итерации, вынося первую итерацию отдельно
            double x = (CountOfIterations > 1) ? h : h/2;
            if (x == h)
            {

            }




            return x;
        }
    }
}
