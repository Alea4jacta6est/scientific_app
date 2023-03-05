FROM python:3.10

ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt

ADD . /scientific_app
WORKDIR /scientific_app
ENV PYTHONPATH /scientific_app

EXPOSE 8501
CMD streamlit run app.py