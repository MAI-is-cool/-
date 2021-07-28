//Class which describe main attributes of residual, such as Name and Value of Residual. Value of residual determines for each step.
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    public class TResidual
    {
        /// <summary>
        /// Название граничной величены
        /// </summary>
        public string Name { get; set; }
        /// <summary>
        /// Значение граничной величены
        /// </summary>
        public double ControlValue { get; set; }
        /// <summary>
        /// Потребная точность невязки
        /// </summary>
        public double Accuracy { get; set; }
        /// <summary>
        /// Текущее значение граничной величены
        /// </summary>
        public double CurrentValue { get; set; }
        /// <summary>
        /// 
        /// </summary>
        public double Residual
        {
            get
            {
                return ControlValue - CurrentValue;
            }
                
        }
    }
}
