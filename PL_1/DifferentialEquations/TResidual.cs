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
        /// Name of residual
        /// </summary>
        public string Name { get; set; }
        /// <summary>
        /// Value of residual
        /// </summary>
        public double Value { get; set; }
    }
}
