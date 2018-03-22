# Jlab Simulation Package for Electron Cooling

## About JSPEC
JSPEC is an open source C++ package for numerical simulations on the electron cooling process, including the intrabeam scattering (IBS) effect, developed at [Jefferson Lab (JLab)](http://www.jlab.org). 



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## How to compile
#### Compile using code blocks IDE

JSPEC is developed using the code blocks IDE. If you are using the same IDE, just open the project file "jspec.cbp" in the cbp folder and build it. Make sure you have a c++ compiler that supports C++11 standard. 

#### Compile using cmake

Users can also use cmake to compile the files (tested in Ubuntu 16.06). In the project folder, run the following commands:

> ` cd build` 
>
> `cmake ..`
>
> `make`



## Acknowledgement

Authors of [**BETACOOL**](http://betacool.jinr.ru/), we learned a lot from BETACOOL. 

Authors of [**muParser**](http://beltoforion.de/article.php?a=muparser),  which is used in building the text-base UI. 

## Contact 

contact the author by hezhang.AT.jlab.org. 