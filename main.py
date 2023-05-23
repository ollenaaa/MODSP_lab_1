import numpy as np
from scipy.integrate import quad
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.io as pio

n = 21
N = 120

# функція f(x)
def f(x):
    return n * np.sin(np.pi * n * x)


# обчислення коефіцієнтів Фур'є
def fourier_coefficient():
    a0 = 2 / np.pi * quad(f, -np.pi/2, np.pi/2)[0]
    ak = []
    bk = []
    for k in range(1, N + 1):
        a = 2 / np.pi * quad(lambda x: f(x) * np.cos(2 * k * x),-np.pi/2, np.pi/2, limit=100)[0]
        b = 2 / np.pi * quad(lambda x: f(x) * np.sin(2 * k * x),0, np.pi, limit=100)[0]
        ak.append(a)
        bk.append(b)
    return a0, ak, bk


# функція для розрахунку наближення функції за допомогою ряду Фур'є
def fourier(x, a0, ak, bk, l):
    return a0/2 + sum([bk[i-1]*np.sin(2*i*x) + ak[i-1]*np.cos(2*i*x) for i in range(l)])


def fourier_approx(x, a0, ak, bk, n):
    s = a0/2
    for k in range(1, n + 1):
        s += ak[k-1]*np.cos(2*k*x) + bk[k-1]*np.sin(2*k*x)
    return s


x = np.linspace(0, np.pi, 2000)
N_values = list(range(1, N+1))
a0, ak, bk = fourier_coefficient()

fig = make_subplots(rows=1, cols=1)

fig.add_trace(go.Scatter(x=x, y=f(x), name="f(x)"), row=1, col=1)
for N in N_values:
    s_values = [fourier_approx(i, a0, ak, bk, N) for i in x]
    label = f"N_{N}(x)"  # підпис для легенди
    fig.add_trace(go.Scatter(x=x, y=s_values, name=label), row=1, col=1)

fig.update_layout(title="Наближення функції за допомогою ряду Фур'є")
fig.show()


fx_fourier = np.array([fourier_approx(i, a0, ak, bk, N) for i in x])
fx = f(x)

ak.insert(0, a0)

error = np.max(np.abs(fx - fx_fourier))
relative_error = error / np.max(np.abs(fx))

trace_fx = go.Scatter(x=x, y=fx, mode='lines', name='функція')
trace_fx_fourier = go.Scatter(x=x, y=fx_fourier, mode='lines', name='наближення')
layout = go.Layout(title=f"Ряд Фур'є (N={N})", xaxis=dict(title='x'), yaxis=dict(title='y'))

fig = go.Figure(data=[trace_fx, trace_fx_fourier], layout=layout)
pio.show(fig)


x_new = np.linspace(-np.pi, 0, 2000)
x_new = np.concatenate((x_new, np.linspace(np.pi, 2*np.pi, 2000)))
fx_fourier = np.array([fourier_approx(i, a0, ak, bk, N) for i in x_new])
fig = go.Figure()
fig.add_trace(go.Scatter(x=x_new, y=fx_fourier, mode='lines'))
fig.show()

i = np.linspace(0, N, N + 1)
ak_plotly = np.array(ak)

fig = go.Figure(data=go.Scatter(x=i, y=ak_plotly, mode='markers+lines'))
fig.update_layout(title='Fourier Coefficients a_k', xaxis_title='k', yaxis_title='a_k')

fig.show()

i = np.linspace(1, N, N)
bk_plotly = np.array(bk)

fig = go.Figure(data=go.Scatter(x=i, y=bk_plotly, mode='markers+lines'))
fig.update_layout(title='Fourier Coefficients b_k', xaxis_title='k', yaxis_title='b_k')

fig.show()

# зберігання у файл
file = open("Lab1.txt", "w")
file.write("a. Порядок " + str(N))
file.write("\nb. Обчислені коефіцієнти тригонометричного ряду Фур'є")
for i in range (len(ak)):
    file.write("\na" + str(i) + " = " + str(ak[i]))
for i in range (len(bk)):
    file.write("\nb" + str(i) + " = " + str(bk[i]))
file.write("\nc. Похибка наближення " + str(relative_error))
file.close()

