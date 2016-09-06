import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in('username', 'api_key')
data = Data([
    Scatter(
        mode='markers',
        stream=Stream(
            maxpoints=50,
            token='your_stream_token_here'
        )
    )
])
plot_url = py.plot(data)