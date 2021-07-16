//Делегат системы дифференциальных уравнений с n переменными и n неизвестными
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Делегат системы дифференциальных уравнений с n переменными и n неизвестными
    /// </summary>
    /// <param name="X">Переменная, по которой идет интегриривание</param>
    /// <param name="Y">Массив неизветных функций, зависящих от x </param>
    /// <returns></returns>
    public delegate double [] AEquation_dN(double X, double[] Y);

}
